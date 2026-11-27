"""Read-only physical-memory process inventory shared by runtime and UI observer."""
import ctypes
import ctypes.util
import os
from pathlib import Path
import sys

class MacProcesses:
    def __init__(self):
        if sys.platform != 'darwin':
            raise RuntimeError('This observer uses macOS physical-footprint accounting')
        self.proc = ctypes.CDLL(ctypes.util.find_library('proc'), use_errno=True)
        self.libc = ctypes.CDLL(None, use_errno=True)
        self.proc.proc_listallpids.argtypes = [ctypes.c_void_p, ctypes.c_int]
        self.proc.proc_listallpids.restype = ctypes.c_int
        self.proc.proc_pid_rusage.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_void_p]
        self.proc.proc_pid_rusage.restype = ctypes.c_int
        self.proc.proc_pidpath.argtypes = [ctypes.c_int, ctypes.c_void_p, ctypes.c_uint32]
        self.proc.proc_pidpath.restype = ctypes.c_int
        self.proc.proc_pidinfo.argtypes = [ctypes.c_int, ctypes.c_int, ctypes.c_uint64, ctypes.c_void_p, ctypes.c_int]
        self.proc.proc_pidinfo.restype = ctypes.c_int
        self.libc.sysctl.argtypes = [ctypes.POINTER(ctypes.c_int), ctypes.c_uint, ctypes.c_void_p,
                                    ctypes.POINTER(ctypes.c_size_t), ctypes.c_void_p, ctypes.c_size_t]
        self.libc.sysctl.restype = ctypes.c_int

    def command(self, pid):
        mib = (ctypes.c_int * 3)(1, 49, pid)  # CTL_KERN, KERN_PROCARGS2, pid
        size = ctypes.c_size_t(262144); raw = ctypes.create_string_buffer(size.value)
        if self.libc.sysctl(mib, 3, raw, ctypes.byref(size), None, 0):
            return None
        argc = ctypes.c_int.from_buffer(raw).value
        values = raw.raw[4:size.value]
        end = values.find(b'\0')
        if end < 0 or not 0 <= argc <= 10000:
            return None
        return [a.decode(errors='replace') for a in values[end+1:].lstrip(b'\0').split(b'\0')[:argc]]

    def ancestry(self, pid):
        rows = []; seen = set()
        for _ in range(16):
            if pid <= 0 or pid in seen:
                break
            seen.add(pid); raw = ctypes.create_string_buffer(256)
            # PROC_PIDTBSDINFO=3; pbi_pid/ppid are uint32 fields 3/4.
            size = self.proc.proc_pidinfo(pid, 3, 0, raw, len(raw))
            if size < 20:
                rows.append(dict(pid=pid, argv=self.command(pid), parent_unreadable=True))
                break
            info = (ctypes.c_uint32 * 5).from_buffer(raw)
            if info[3] != pid:
                raise RuntimeError('BSD process identity did not match observed PID')
            parent = int(info[4]); rows.append(dict(pid=pid, parent_pid=parent, argv=self.command(pid)))
            pid = parent
        return rows

    def sample(self):
        count = self.proc.proc_listallpids(None, 0)
        if count <= 0:
            raise RuntimeError('Cannot enumerate processes')
        pids = (ctypes.c_int * (count + 1024))()
        count = self.proc.proc_listallpids(pids, ctypes.sizeof(pids))
        if count <= 0:
            raise RuntimeError('Cannot read process list')
        rows, unreadable = [], []
        for pid in pids[:count]:
            if pid <= 0:
                continue
            raw = ctypes.create_string_buffer(96)  # rusage_info_v0
            if self.proc.proc_pid_rusage(pid, 0, raw):
                unreadable.append(dict(pid=pid, errno=ctypes.get_errno()))
                continue
            values = (ctypes.c_uint64 * 10).from_buffer(raw, 16)
            path = ctypes.create_string_buffer(4096)
            self.proc.proc_pidpath(pid, path, len(path))
            rows.append(dict(pid=pid, process_start_abstime=int(values[8]),
                             rss_bytes=int(values[6]), physical_footprint_bytes=int(values[7]),
                             executable=path.value.decode(errors='replace')))
        return rows, unreadable



class LinuxProcesses:
    """Linux exposes RSS + VmSwap rather than macOS physical footprint."""
    def command(self, pid):
        try:return [x.decode(errors="replace") for x in (Path("/proc")/str(pid)/"cmdline").read_bytes().split(b"\0") if x]
        except OSError:return None

    def ancestry(self, pid):
        rows=[];seen=set()
        for _ in range(16):
            if pid<=0 or pid in seen:break
            seen.add(pid)
            try:
                text=(Path("/proc")/str(pid)/"stat").read_text();fields=text[text.rfind(")")+2:].split();parent=int(fields[1])
            except OSError:
                rows.append(dict(pid=pid,argv=self.command(pid),parent_unreadable=True));break
            rows.append(dict(pid=pid,parent_pid=parent,argv=self.command(pid)));pid=parent
        return rows

    def sample(self):
        rows=[];unreadable=[];page=os.sysconf("SC_PAGE_SIZE")
        for directory in Path("/proc").iterdir():
            if not directory.name.isdigit():continue
            pid=int(directory.name)
            try:
                text=(directory/"stat").read_text();fields=text[text.rfind(")")+2:].split();rss=int(fields[21])*page
                status=(directory/"status").read_text().splitlines();swap=[line for line in status if line.startswith("VmSwap:")];swapped=int(swap[0].split()[1])*1024 if swap else 0
                try:executable=os.readlink(directory/"exe")
                except OSError:executable=""
                rows.append(dict(pid=pid,process_start_abstime=int(fields[19]),rss_bytes=rss,physical_footprint_bytes=rss+swapped,executable=executable))
            except OSError as error:unreadable.append(dict(pid=pid,errno=error.errno))
        return rows,unreadable


def process_observer():
    if sys.platform=="darwin":return MacProcesses()
    if sys.platform.startswith("linux"):return LinuxProcesses()
    raise RuntimeError("All-process memory observation supports macOS and Linux only")
