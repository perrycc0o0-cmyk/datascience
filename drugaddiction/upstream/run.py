#!/usr/bin/env python3
"""Run the study-specific snRNA-seq upstream stages with local paths."""

import argparse
import json
from pathlib import Path
import shutil
import subprocess


ROOT = Path(__file__).resolve().parent


def config(name):
    return json.loads((ROOT / 'config' / name).read_text(encoding='utf-8'))


def put(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        if path.read_bytes() != data:
            raise SystemExit(f'Existing run input differs: {path}; choose a new run directory.')
    else:
        path.write_bytes(data)


def write_config(root, name, value):
    put(root / 'config' / name, (json.dumps(value, indent=2) + '\n').encode('utf-8'))


def copy_scripts(root, stage):
    for path in sorted((ROOT / stage).glob('*.R')):
        put(root / 'scripts' / path.name, path.read_bytes())


def rscript(exe, root, name, *args):
    path = root / 'scripts' / name
    print(f'Running {name}', flush=True)
    subprocess.run([exe, '--vanilla', str(path), *map(str, args)], check=True)


def require(path):
    if not path.is_file():
        raise SystemExit(f'Required input not found: {path}')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('step', choices=('preprocess', 'global', 'neurons', 'prepare'))
    parser.add_argument('--run-dir', required=True, type=Path, help='Private work directory for all stage outputs')
    parser.add_argument('--input-dir', type=Path, help='Parent of nine sample matrix directories (preprocess only)')
    parser.add_argument('--rscript', default='Rscript', help='Rscript executable or executable name on PATH')
    parser.add_argument('--fontconfig-file', type=Path, help='Optional local Fontconfig configuration for Arial')
    parser.add_argument('--resume', action='store_true', help='Resume the per-library checkpointed preprocess step')
    a = parser.parse_args()
    if a.step != 'preprocess' and (a.input_dir or a.resume):
        parser.error('--input-dir and --resume only apply to preprocess')
    if a.step not in ('preprocess', 'neurons') and a.fontconfig_file:
        parser.error('--fontconfig-file only applies when creating the preprocess or neurons configuration')
    if a.step == 'preprocess' and a.input_dir is None:
        parser.error('preprocess requires --input-dir')
    if a.fontconfig_file:
        require(a.fontconfig_file.resolve())
    executable = shutil.which(a.rscript)
    if not executable:
        parser.error(f'Rscript executable not found: {a.rscript}')
    root = a.run_dir.resolve()
    if root == ROOT.parent or ROOT.parent in root.parents:
        parser.error('Use a work directory outside the public code folder.')
    stage1, stage2, stage3 = (root / f'stage{i}' for i in range(1, 4))

    if a.step == 'preprocess':
        cfg = config('processing_config.json')
        cfg['input_root'] = str(a.input_dir.resolve())
        cfg['runtime']['fontconfig_file'] = str(a.fontconfig_file.resolve()) if a.fontconfig_file else None
        write_config(stage1, 'processing_config.json', cfg)
        put(stage1 / 'config' / 'sample_naming.json', (ROOT / 'config' / 'sample_naming.json').read_bytes())
        copy_scripts(stage1, 'stage1')
        rscript(executable, stage1, '00_runtime_preflight.R', stage1)
        rscript(executable, stage1, '01_process_libraries.R', '--root', stage1, *(['--resume'] if a.resume else []))
    elif a.step == 'global':
        require(stage1 / '02_qc_doublet' / 'attrition_summary.csv')
        for name in ('02_global_celltype.R', '03_plot_global_review.R',
                     '06_match_six_major_classes.R', '07_plot_six_class_review.R'):
            rscript(executable, stage1, name, '--root', stage1)
    elif a.step == 'neurons':
        source = stage1 / '05_six_class_review' / 'global_celltype_six_class.rds'
        require(source)
        cfg = config('stage2_config.json')
        cfg['parent_stage1_root'] = str(stage1)
        cfg['input_rds'] = str(source)
        cfg['runtime']['fontconfig_file'] = str(a.fontconfig_file.resolve()) if a.fontconfig_file else None
        write_config(stage2, 'stage2_config.json', cfg)
        put(stage2 / 'config' / 'subtype_markers_and_style.json', (ROOT / 'config' / 'subtype_markers_and_style.json').read_bytes())
        copy_scripts(stage2, 'stage2')
        for name in ('00_runtime_preflight.R', '01_recluster_neurons.R', '02_annotate_neuron_subtypes.R'):
            rscript(executable, stage2, name, '--root', stage2)
        rscript(executable, stage2, '03_plot_neuron_review.R', '--root', stage2,
                '--runtime-config', stage2 / 'config' / 'stage2_config.json')
    else:
        require(stage2 / '03_annotation' / 'Exc_cell_annotations_UMAP.csv.gz')
        require(stage2 / '03_annotation' / 'Inh_cell_annotations_UMAP.csv.gz')
        if (stage3 / '01_inputs').exists():
            raise SystemExit(f'Input preparation has already started in {stage3}; preserve it before starting a new preparation.')
        cfg = config('analysis_contract.json')
        cfg['base'] = str(root)
        write_config(stage3, 'analysis_contract.json', cfg)
        copy_scripts(stage3, 'stage3')
        rscript(executable, stage3, '00_apply_accepted_subtype_override.R', stage2)
        rscript(executable, stage3, '01_prepare_analysis_inputs.R', stage3)


if __name__ == '__main__':
    main()
