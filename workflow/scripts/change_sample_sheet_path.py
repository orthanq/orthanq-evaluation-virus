import yaml
import os

with open(snakemake.log[0], "w") as f:
    sys.stderr = sys.stdout = f

    class MyDumper(yaml.Dumper):

        def increase_indent(self, flow=False, indentless=False):
            return super(MyDumper, self).increase_indent(flow, False)
    
    #update pep config to have the poath for newly generated sample sheet
    with open(snakemake.input.pepfile) as f:
        doc = yaml.safe_load(f)

        doc['sample_table'] = os.path.basename(snakemake.params.sample_sheet_path)
        
    with open(snakemake.output.pep_config, 'w') as f:
        yaml.dump(doc, f, sort_keys=False, Dumper=MyDumper, default_flow_style=False)       


    #update main config to have the right file path for pep file
    with open(snakemake.input.config) as f:
        doc = yaml.safe_load(f)

        doc['pepfile'] = "config/pep/" + os.path.basename(snakemake.params.pepfile_path)
        print(doc)
    with open(snakemake.output.main_config, 'w') as f:
        yaml.dump(doc, f, sort_keys=False, Dumper=MyDumper, default_flow_style=False)
