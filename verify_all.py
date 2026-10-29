"""Verify all manuscript certificate claims with Python 3.10+ standard library."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib, json, subprocess, sys, time
ROOT=Path(__file__).resolve().parent
def main():
    started=time.monotonic();out=ROOT/'verification';out.mkdir(exist_ok=True)
    checked=0
    for manifest in sorted((ROOT/'certificates').glob('*/MANIFEST.sha256')):
        for line in manifest.read_text().splitlines():
            expected,rel=line.split('  ',1);p=manifest.parent/rel
            if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=expected:
                raise ValueError('Frozen-file mismatch: '+str(p.relative_to(ROOT)))
            checked+=1
    jobs=[
      ('frontier','illumination-testing-frontier','verify_refinement.py',['--check-faults']),
      ('stability','illumination-testing-stability','verify_stability.py',[]),
      ('stability-faults','illumination-testing-stability','check_faults.py',[]),
      ('prediction','illumination-testing-prediction','verify_prediction.py',['--check-faults'])]
    reports={}
    for name,folder,script,extra in jobs:
        report=out/(name+'.json');log=out/(name+'.log')
        print('Checking '+name+' ...',flush=True)
        with log.open('w',encoding='utf-8') as f:
            subprocess.run([sys.executable,'-I','-O','-B',str(ROOT/'certificates'/folder/script),*extra,'--output',str(report)],cwd=ROOT,stdout=f,stderr=subprocess.STDOUT,check=True)
        r=json.loads(report.read_text());
        if r.get('status')!='PASS':raise ValueError('Failed report: '+name)
        reports[name]={'status':'PASS','report':report.relative_to(ROOT).as_posix()}
        print(name+': PASS',flush=True)
    result=dict(status='PASS',createdAtUtc=datetime.now(timezone.utc).isoformat(),
                checkedFrozenFiles=checked,optimizerCalls=0,reports=reports,
                elapsedSeconds=round(time.monotonic()-started,3))
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
if __name__=='__main__':main()
