"""Reject invalid stability claims with exact reconstructed certificates."""
from pathlib import Path
from fractions import Fraction as F
import argparse, copy, importlib.util, json
HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('stability',HERE/'verify_stability.py')
s=importlib.util.module_from_spec(spec);spec.loader.exec_module(s)

def main():
    p=s.v.read(HERE/'proof.json');data=s.prepare(p);s.check_claims(p,data);faults=[]
    def reject(name,operation,message):
        try:operation()
        except s.v.InvalidCertificate as error:
            s.require(message in str(error),'Unexpected rejection reason')
            faults.append(dict(name=name,status='REJECTED_AS_EXPECTED',reason=str(error)))
        else:raise RuntimeError('Invalid claim was accepted: '+name)
    q=copy.deepcopy(p);q['coreProofSHA256']='0'*64
    reject('changed core binding',lambda:s.check_configuration(q),'Core proof binding')
    q=copy.deepcopy(p);q['certificates'][0]['Uupper']='0'
    reject('zero model residual norm',lambda:s.check_claims(q,data),'Residual norm understated')
    q=copy.deepcopy(p);q['certificates'][0]['Supper']='0'
    reject('zero shared block norm',lambda:s.check_claims(q,data),'Shared block norm understated')
    q=copy.deepcopy(p);q['rhoAxis']=str(2*F(p['rhoAxis']))
    reject('doubled model allowance',lambda:s.check_claims(q,data),'Norm triangle unsafe')
    q=copy.deepcopy(p);q['kappaAxis']=str(2*F(p['kappaAxis']))
    reject('doubled shape allowance',lambda:s.check_claims(q,data),'Shape triangle unsafe')
    reject('combining both full norm-axis allowances',lambda:s.certify_norm_point(data,F(p['muAxis']),F(p['rhoAxis'])),'Point outside norm guarantee')
    q=copy.deepcopy(p);q['assumptions']['normalization']='nominal-white-Y'
    reject('wrong normalization assumption',lambda:s.check_configuration(q),'Stability assumptions changed')
    return dict(status='PASS',optimizerCalls=0,analyticControls=s.controls(),faultCases=faults)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--output',type=Path);args=parser.parse_args()
    result=main()
    if args.output:args.output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))
