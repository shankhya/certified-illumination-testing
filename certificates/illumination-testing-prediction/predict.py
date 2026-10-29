"""Exact evaluation of the certified affine coefficients, using rational input."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json

def rationals(values,length):
    if len(values)!=length or any(not isinstance(x,(str,int)) or isinstance(x,bool) for x in values):
        raise ValueError('Use exact rational/decimal strings or integers with the declared dimensions')
    return [F(x) for x in values]

def predict(model,observation):
    y=rationals(observation['y'],15);wM=rationals(observation['testWhite'],15);wQ=rationals(observation['targetWhite'],9)
    V=[rationals(row,15) for row in model['V']]
    if len(V)!=9:raise ValueError('Nine target coordinates required')
    return [sum((a*(t-w/2) for a,t,w in zip(row,y,wM)),F(0))+q/2 for row,q in zip(V,wQ)]

def main():
    parser=argparse.ArgumentParser();parser.add_argument('observation');parser.add_argument('--output');args=parser.parse_args()
    model=json.loads((Path(__file__).resolve().parent/'proof.json').read_text())
    observed=json.loads(Path(args.observation).read_text());value=predict(model,observed)
    out=dict(targetResponseExact=list(map(str,value)),targetResponseApproximate=list(map(float,value)),
             note='Exact algebraic evaluation. The error guarantee additionally requires the verified class, noise, kernel and white-calibration assumptions.')
    text=json.dumps(out,indent=2)+'\n'
    if args.output:Path(args.output).write_text(text)
    else:print(text)
if __name__=='__main__':main()
