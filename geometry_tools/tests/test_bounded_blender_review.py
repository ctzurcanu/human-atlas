"""Exercise queue selection/lifetime without importing or launching Blender."""
import ast
import hashlib
import json
import os
from pathlib import Path
import tempfile
import time
import unittest


class BoundedBlenderReviewTests(unittest.TestCase):
    def environment(self,root,request):
        root=root.resolve()
        source=ast.parse(Path('geometry_tools/blender_bridge.py').read_text())
        names={'_ga_step','_ga_atomic','_ga_path','_ga_sha'}
        body=[n for n in source.body if isinstance(n,ast.FunctionDef)and n.name in names]
        env=dict(Path=Path,json=json,hashlib=hashlib,os=os,time=time,_ga_queue=root,_ga_request_id=request,_ga_bridge_deadline=time.monotonic()+60)
        exec(compile(ast.Module(body=body,type_ignores=[]),'bounded-review-test','exec'),env)
        calls=[];env['_ga_render']=lambda task:calls.append(task['id'])or{'request_id':task['id'],'state':'rendered'}
        return env,calls

    def queue(self,root,request):
        root=root.resolve()
        folder=root/request;folder.mkdir();p=folder/'request.json';p.write_text(json.dumps({'id':request,'action':'render','output':str(folder)}))
        pending=root/'pending';pending.mkdir(exist_ok=True)
        (pending/(request+'.json')).write_text(json.dumps({'request':str(p),'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}))

    def test_only_selected_request_runs_and_timer_stops(self):
        with tempfile.TemporaryDirectory()as temp:
            root=Path(temp);first='a'*64;other='b'*64;self.queue(root,first);self.queue(root,other)
            env,calls=self.environment(root,first)
            self.assertIsNone(env['_ga_step']());self.assertEqual(calls,[first])
            self.assertTrue((root/'pending'/(other+'.json')).exists())
            self.assertEqual(json.loads((root/'bridge.json').read_text())['state'],'stopped')

    def test_expired_session_leaves_pending_work_untouched(self):
        with tempfile.TemporaryDirectory()as temp:
            root=Path(temp);request='a'*64;self.queue(root,request);env,calls=self.environment(root,request);env['_ga_bridge_deadline']=time.monotonic()-1
            self.assertIsNone(env['_ga_step']());self.assertEqual(calls,[]);self.assertTrue((root/'pending'/(request+'.json')).exists())

    def test_missing_explicit_identity_never_replays_old_queue(self):
        with tempfile.TemporaryDirectory()as temp:
            root=Path(temp);self.queue(root,'a'*64);env,calls=self.environment(root,None)
            self.assertIsNone(env['_ga_step']());self.assertEqual(calls,[])
