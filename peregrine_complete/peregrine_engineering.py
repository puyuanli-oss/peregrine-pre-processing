#!/usr/bin/env python3
"""Auditable PEREGRINE orchestrator. Does not modify legacy linking algorithms."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
OUTPUTS = {'eQTL':'eQTL/linksDBnumeqtl', 'hierarchicalTAD':'heirarchicalTAD/PSYCHIClinksDB', 'ChIA-PET':'chia_pet/linksDBchia', 'TAD':'tad/linksDBtad'}

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda:f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()

def validate_links(path, assay=None):
    """Check headerless link contract; empty results are valid but missing are not."""
    if not path.is_file(): raise FileNotFoundError(f'Missing result: {path}')
    count=0; examples=[]
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            count+=1
            fields=line.rstrip('\r\n').split('\t')
            if len(fields)!=4 or any(not x for x in fields) or (assay is not None and fields[3]!=assay):
                if len(examples)<10: examples.append(f'{count}: invalid 4-field record or assay: {line[:120].rstrip()}')
    if examples: raise ValueError(f'{path}: '+ '; '.join(examples))
    return {'path':str(path),'rows':count,'bytes':path.stat().st_size,'sha256':digest(path)}

def validate_manifest(path):
    if not path.is_file(): raise FileNotFoundError(f'Missing ChIA manifest: {path}')
    with path.open(newline='',encoding='utf-8') as f:
        reader=csv.DictReader(f,delimiter='\t')
        expected={'sample_id','cell_type','target_protein','mango_file'}
        if not expected.issubset(reader.fieldnames or []): raise ValueError('ChIA manifest missing required columns')
        seen=set(); samples=[]
        for row in reader:
            sid=row['sample_id']
            if not sid or sid in seen: raise ValueError('Duplicate or blank ChIA sample ID')
            seen.add(sid)
            path_value=Path(row['mango_file'])
            sample=path_value if path_value.is_absolute() else path.parent/path_value
            if not sample.is_file() or sample.stat().st_size==0: raise FileNotFoundError(f'ChIA sample missing or empty: {sample}')
            samples.append({'sample':sid,'file':str(sample.resolve()),'sha256':digest(sample)})
    if not samples: raise ValueError('Empty ChIA manifest')
    return samples

def run(command, cwd, log):
    log.parent.mkdir(parents=True,exist_ok=True)
    with log.open('w',encoding='utf-8') as out:
        out.write('Command: '+repr([str(x) for x in command])+'\n'); out.flush()
        result=subprocess.run([str(x) for x in command],cwd=cwd,stdout=out,stderr=subprocess.STDOUT,check=False)
    if result.returncode: raise RuntimeError(f'Command failed ({result.returncode}): {command}; inspect {log}')

def main(argv=None):
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('command',choices=['preflight','preprocess','run','validate'])
    p.add_argument('--output',type=Path,default=Path('results'))
    p.add_argument('--manifest',type=Path,default=ROOT/'chia_pet/samples.tsv')
    p.add_argument('--genes',type=Path,default=ROOT/'chia_pet/hg19genes.txt')
    p.add_argument('--gtf',type=Path)
    p.add_argument('--offline',action='store_true')
    p.add_argument('--skip-preprocess',action='store_true')
    p.add_argument('--allow-missing-chia',action='store_true',help='Explicitly disable ChIA when its raw files are unavailable')
    args=p.parse_args(argv)
    output=args.output.resolve(); manifest=args.manifest.resolve(); genes=args.genes.resolve()
    audit=output/'engineering'; audit.mkdir(parents=True,exist_ok=True)
    report={'schema_version':1,'started_utc':datetime.now(timezone.utc).isoformat(),'command':args.command,'stages':[],'status':'running'}
    def stage(name,cmd):
        logfile=audit/(f'{len(report["stages"])+1:02d}_{name}.log')
        item={'name':name,'command':[str(x) for x in cmd],'log':str(logfile),'status':'running'}
        report['stages'].append(item)
        try: run(cmd,ROOT,logfile); item['status']='passed'
        except Exception: item['status']='failed'; raise
    try:
        if args.command=='validate':
            report['outputs']={name:validate_links(output/path,'1' if name=='ChIA-PET' else None) for name,path in OUTPUTS.items()}
        else:
            if not (ROOT/'CREbedDBenhancers_10092018').is_file():
                raise FileNotFoundError('Required curated CREbedDBenhancers_10092018 is not included; supply the original file (do not fabricate it).')
            if not genes.is_file(): raise FileNotFoundError(f'ChIA gene table missing: {genes}')
            chia=True
            try: report['chia_inputs']=validate_manifest(manifest)
            except (FileNotFoundError,ValueError):
                if not args.allow_missing_chia: raise
                chia=False; report['chia_skipped']='Explicitly allowed missing ChIA input; not a complete four-assay run'
            if args.command=='preflight':
                if shutil.which('bedtools') is None: raise RuntimeError('bedtools missing from PATH')
                stage('preflight',[sys.executable,'verify_peregrine.py','--repo',ROOT,'--scope','chia' if chia else 'final','--output',output,'--check-inputs'] if chia else [sys.executable,'-m','unittest','tests_preprocessing'])
            else:
                if not args.skip_preprocess:
                    cmd=[sys.executable,'-m','preprocess.cli','all']
                    if args.offline: cmd.append('--offline')
                    if args.gtf: cmd+=['--gtf',args.gtf]
                    stage('preprocessing',cmd)
                if args.command=='run':
                    if not chia: raise RuntimeError('Four-assay full run requires ChIA inputs; cannot silently skip ChIA')
                    stage('linking',['bash','run_peregrine.sh',output,'--chia-manifest',manifest,'--chia-genes',genes])
                    report['outputs']={name:validate_links(output/path,'1' if name=='ChIA-PET' else None) for name,path in OUTPUTS.items()}
        report['status']='passed'
    except Exception as exc:
        report['status']='failed'; report['error']=f'{type(exc).__name__}: {exc}'
        print(report['error'],file=sys.stderr)
    finally:
        report['finished_utc']=datetime.now(timezone.utc).isoformat()
        (audit/'report.json').write_text(json.dumps(report,indent=2)+'\n')
        print('Audit report:',audit/'report.json')
    return 0 if report['status']=='passed' else 1

if __name__=='__main__': sys.exit(main())
