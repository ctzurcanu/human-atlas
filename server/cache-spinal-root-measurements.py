"""Cache primary morphometric data for procedural rootlet/sleeve planning.

One published female case is reference evidence, not universal anatomy or
donor-specific truth for the atlas. No third-party geometry is downloaded.
"""
from pathlib import Path
import argparse,hashlib,json,time,urllib.request,xml.etree.ElementTree as ET

ap=argparse.ArgumentParser();ap.add_argument('--output',required=True);args=ap.parse_args();start=time.monotonic();out=Path(args.output).resolve();out.mkdir(parents=True,exist_ok=False)
url='https://www.ebi.ac.uk/europepmc/webservices/rest/PMC5738087/fullTextXML'
request=urllib.request.Request(url,headers={'User-Agent':'Human-Atlas-local-reference-review','Accept':'application/xml'})
with urllib.request.urlopen(request,timeout=10)as response:data=response.read(2*1024*1024+1)
assert len(data)<=2*1024*1024;root=ET.fromstring(data)
assert any(e.text=='10.1186/s12987-017-0085-y'for e in root.iter('article-id'))
tables=[t for t in root.iter('table-wrap')if ''.join(t.find('label').itertext()).strip()in['Table 1','Table 1']];assert len(tables)==1
rows=[]
for tr in tables[0].findall('.//tbody/tr'):
    cells=[''.join(e.itertext()).strip()for e in list(tr)if e.tag in['td','th']]
    if len(cells)!=10:continue
    name=cells[0];numbers=[float(''.join(s.split()).replace('−','-'))for s in cells[1:]]
    def measurement(values):
        angle,line,diameter,count=values;assert count==int(count)and count>0
        return dict(descending_angle_deg=angle,radicular_line_mm=line,diameter_mm=diameter,rootlet_count_per_side=int(count))
    rows.append(dict(level=name,dural_exit_distance_from_foramen_magnum_mm=numbers[0],dorsal=measurement(numbers[1:5]),ventral=measurement(numbers[5:9])))
expected=[*('C'+str(i)for i in range(1,9)),*('T'+str(i)for i in range(1,13)),*('L'+str(i)for i in range(1,6)),*('S'+str(i)for i in range(1,6))]
assert [r['level']for r in rows]==expected
result=dict(doi='10.1186/s12987-017-0085-y',primary_article='https://link.springer.com/article/10.1186/s12987-017-0085-y',source_url=url,
    source_xml_sha256=hashlib.sha256(data).hexdigest(),reference='Sass et al. Fluids and Barriers of the CNS 14, 36 (2017), Table 1',
    source_scope='Single adult female case; rootlets semi-idealized from MRI and published cadaveric measurements; values are not universal or atlas-donor-specific.',
    model_geometry_license='CC BY-SA 4.0; no geometry downloaded or installed',article_license='CC BY 4.0; publisher applies CC0 to article data unless otherwise stated',
    measurements=rows,coccygeal_limit='Co1 is not tabulated; the study describes a continuous coccygeal stem before the intrathecal exit. Do not invent a tabulated Co1 dorsal/ventral pair.',
    anatomical_acceptance=False,installed=False,elapsed_seconds=time.monotonic()-start)
(out/'root-measurements.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(dict(levels=len(rows),source_xml_sha256=result['source_xml_sha256'],geometry_downloaded=False,elapsed_seconds=result['elapsed_seconds'])),flush=True)
