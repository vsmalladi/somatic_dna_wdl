import sys
import json
import re
import logging as log
import pandas as pd
import numpy as np
import Colorer
#plots
import plotly.express as px
import plotly.offline as po
# import colorlover as cl
from plotly.subplots import SubplotRef
from plotly import tools
from plotly.subplots import make_subplots
from collections import OrderedDict
#custom
import compose
import ploter


pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class PlotRuntime():
    '''Plot results from outputMetrics.csv'''
    def __init__(self, 
                 output_info_file,
                 metrics_file,
                 plot_file=False):
        self.set_colors()
        if not plot_file:
            self.out_md = metrics_file.replace('.csv', '.md')
            self.header_file = metrics_file.replace('.csv', '.header.txt')
            self.plot_file = metrics_file.replace('.csv', '.html')
        else:
            self.out_md = plot_file.replace('.html', '.md')
            self.header_file = plot_file.replace('.html', '.header.txt')
            self.plot_file = plot_file
        self.output_info = self.load_input(output_info_file)
        self.metadata = self.load_metadata(metrics_file)
        # intro figs
        self.summary_table = self.get_table()
        self.gather_summary()
        self.plot_summary()
        # main figs
        self.fig = self.plot_by_steps(data_steps=self.metadata,
                                      x="workflow_name",
                                      ys=["sample_task_core_h", 
                                          "max_mem_g", 
                                          "sample_task_run_time_h"],
                                      color="task_call_name",
                                      hover_data=['id'],
                                      color_discrete_sequence=self.colors_set2,
                                      title='Tasks: resource usage summary plot',
                                      labels={'sample_task_run_time_h' : 'Max Task runtime(h)',
                                              'max_mem_g' : 'Mem (G)',
                                              'sample_task_core_h' : 'Task core hours',
                                              'workflow_name' : 'Sub-workflow',
                                              'task_call_name' : 'Task'},
                                      x_title="Pipeline tasks")
        self.fig2 = self.plot_by_steps(data_steps=self.metadata,
                                      x="workflow_name",
                                      ys=["sample_subworkflow_core_h", 
                                          "subworkflow_max_mem_g", 
                                          "sample_subworkflow_run_time_h"],
                                      color="workflow_name",
                                      hover_data=['id'],
                                      color_discrete_sequence=self.colors_set1,
                                      labels={'sample_subworkflow_run_time_h' : 'Max sub-workflow runtime(h)',
                                              'subworkflow_max_mem_g' : 'Mem (G)',
                                              'sample_subworkflow_core_h' : 'Sub-workflow core hours',
                                              'workflow_name' : 'Sub-workflow',
                                              'task_call_name' : 'Task'},
                                      title='Subworkflows: resource usage summary plot',
                                      x_title="Pipeline subworkflows")
        self.fig3 = self.plot_by_steps(data_steps=self.metadata,
                                      x="id",
                                      ys=["sample_workflow_core_h", 
                                          "workflow_max_mem_g", 
                                          "sample_workflow_run_time_h"],
                                      color="id",
                                      hover_data=['id'],
                                      color_discrete_sequence=self.colors,
                                      labels={'sample_workflow_run_time_h' : 'Max workflow runtime(h)',
                                              'workflow_max_mem_g' : 'Mem (G)',
                                              'sample_workflow_core_h' : 'Workflow core hours',
                                              'workflow_name' : 'Sub-workflow',
                                              'task_call_name' : 'Task'},
                                      title='Workflows: resource usage summary plot',
                                      x_title="Full pipeline")
        
    def set_colors(self):
        self.colors = px.colors.qualitative.Dark24
        self.map = {'prep_flowcell' : self.colors[0], 'merge' : self.colors[1]}
        self.colors_set1 = px.colors.qualitative.Set1
        self.colors_set2 = px.colors.qualitative.Alphabet
            
    def load_input(self, output_info_file):
        with open(output_info_file) as output_info_object:
                output_info = json.load(output_info_object)
        return output_info
    
    def load_metadata(self, metrics_file):
        metadata = pd.read_csv(metrics_file)
        return metadata
    
    def gather_summary(self):
        '''Add details about percent of instances that are preemptible'''
        print(self.metadata.shape, self.metadata.drop_duplicates('instance_name').shape)
        # exit status of non-zero exits
        exit_status = self.metadata[['task_call_name', 
                                     'execution_status']].value_counts().to_frame().reset_index()
        self.exit_status = exit_status[exit_status.execution_status != 'Done']
        self.exit_status.columns = ['Task', 'Exit status', 'Count']
        # preemptible percentage
        self.preemptible_percent = self.metadata[self.metadata.preemptible].shape[0] / float(self.metadata.shape[0])
        # disk type
        self.metadata['disk_type']= self.metadata.apply(lambda row: row.disk_types.replace('[\'', '').replace('\']', ''), axis=1)
        self.disk_type = self.metadata.groupby(['id', 'disk_type']).sample_task_run_time_h.sum().reset_index()
        self.disk_type.columns = ['Id', 'Disk type', 'Wall clock (h)']
        
    def plot_summary(self):
        fig = px.box(self.disk_type, x='Disk type', y='Wall clock (h)', points="all",
                          color_discrete_sequence=self.colors_set2,
                          color='Id')
        for trace in fig.data:
            trace['showlegend'] = True
            trace['pointpos'] = 0
        fig.update_layout(xaxis_type='category',
                          title_text='Runtime by disk type')
        fig.update_xaxes(title_text='')
        fig.update_yaxes(title_text="Wall clock (h)")
        self.disk_type_plot = ploter.Fig(fig)
        self.exit_status_table = ploter.make_html_table(self.exit_status,
                                                        name='exit_status',
                                                        col_widths=['50%', '30%', '20%'],
                                                        searching=True,
                                                        scroll=True)
        
    def write_summary(self):
        '''Write about details including percent of instances that are preemptible'''
        lines = []
        lines += [str(int(self.preemptible_percent)) + '% of instances were preemptible.',
                  '']
        
    
    def get_table(self):
        '''Make workflow resource usage summary table'''
        table_data = self.metadata[['id', 'sample_workflow_run_time_h']].copy()
        table_data = table_data.sort_values(['sample_workflow_run_time_h']).copy()
        table_data.columns = ['Id', 'Wallclock (h)']
        summary_table = ploter.make_html_table(table_data.drop_duplicates(),
                                               name='runtime_overview',
                                               col_widths=['80%','20%'],
                                               searching=True,
                                               scroll=True)
        return summary_table
    
    def plot_by_steps(self, 
                      data_steps, 
                      x="workflow_name",
                      ys=["sample_subworkflow_core_h", "mem_total_gb", "sample_task_run_time_h"],
                      color="task_call_name",
                      hover_data=['id'],
                      color_discrete_sequence=False, 
                      color_discrete_map=False, 
                      title='',
                      labels={},
                      x_title="Pipeline subworkflow"):
        ''' Plot Task CPU, Mem and Wallclock'''
        category_orders = {}
        data_steps = data_steps.drop_duplicates(subset=[x, color] + ys).copy()
        for y in ys:
            grouped = data_steps.groupby(x)
            category_orders[y] = grouped[y].agg(['max']).reset_index().sort_values('max')[x].unique().tolist()
        if color_discrete_sequence:
            fig1 = px.box(data_steps, x=x, y=ys[0], points="all",
                          hover_data=hover_data,
                          color_discrete_sequence=color_discrete_sequence,
                          labels=labels,
                          color=color, category_orders=category_orders)
            fig2 = px.box(data_steps, x=x, y=ys[1], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
            fig3 = px.box(data_steps, x=x, y=ys[2], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
        else:
            fig1 = px.box(data_steps, x=x, y=ys[0], points="all",
                          color_discrete_map=color_discrete_map,
                          hover_data=hover_data,
                          labels=labels,
                          color=color, category_orders=category_orders)
            fig2 = px.box(data_steps, x=x, y=ys[1], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_map=color_discrete_map,
                          color=color, category_orders=category_orders)
            fig3 = px.box(data_steps, x=x, y=ys[2], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_map=color_discrete_map,
                          color=color, category_orders=category_orders)
        fig = make_subplots(rows=1,
                            cols=3,
                            subplot_titles=("CPU", "Mem", "Wall clock"))
        for trace in fig1.data:
            trace['showlegend'] = False
            trace['pointpos'] = 0
            fig.add_trace(
                          trace,
                          row=1, col=1
                         )
        for trace in fig2.data:
            trace['showlegend'] = False
            trace['pointpos'] = 0
            fig.add_trace(
                          trace,
                          row=1, col=2
                       )
        for trace in fig3.data:
            trace['showlegend'] = True
            trace['pointpos'] = 0
            fig.add_trace(
                          trace,
                          row=1, col=3
                       )
        fig.update_layout(xaxis_type='category',
                          title_text=title)
        fig.update_xaxes(title_text='')
        fig.update_yaxes(title_text="CPU Time (h)", row=1, col=1)
        fig.update_yaxes(title_text="Mem (G)", row=1, col=2)
        fig.update_yaxes(title_text="Wall clock (h)", row=1, col=3)
        fig.update_xaxes(title_text="",
                        categoryarray=category_orders[ys[0]],
                        categoryorder='array',
                        row=1, col=1)
        fig.update_xaxes(title_text=x_title,
                        categoryarray=category_orders[ys[1]],
                        categoryorder='array',
                         row=1, col=2)
        fig.update_xaxes(title_text="",
                        categoryarray=category_orders[ys[2]],
                        categoryorder='array',
                         row=1, col=3)
        return ploter.Fig(fig)
        
    
def make_files(results, appendix=False):
    '''
        Print out plots made for report.
    '''
    with open(results.out_md, 'w') as report:
        header_file = results.header_file
        #  =======================
        #  Prep Content
        #  =======================
        contents = compose.Content(tumor_only=False,
                                   library='WGS',
                                   internal=False)
        #  =======================
        #  Section Logo and basic info
        #  =======================
        basic_info = ''.join(['\n\nRuntime metrics: v7 pipeline',
                              '\n\nAverage runtime for different subworkflows of the v7 pipeline are shown below..\n\n'  
#                               Pre-processing includes "prep_flowcell" (alignment, short alignment filtering and fixmate), "merge" (duplicate marking and base quality score recalibration). Calling includes somatic and germline variant calling, merging, filtering, annotation, HLA typing, MSI classification and signature estimation.\n\n'
                              ])
#        basic_info = ''.join(['\n\n<b>Tumor :</b> ' + results.tumor +
#                              '\n<b>Normal :</b> ' + results.normal +
#                              '\n\n<b>Project :</b> ' + results.project_name])
#        if results.organism in ['Mouse']:
#            basic_info += '\n\n<b>Organism :</b> ' + results.organism
        logo = '<img src="NY-Genome-Logo.jpg" class="center">' + ('\n' * 4)
        section_content = ''.join([logo, basic_info, '\n\n'])
        all_section_content = contents.add_line(content=section_content,
                                               libraries=['WGS', 'Exome'],
                                               types=['tumor_only', 'paired'],
                                               audiences=['internal', 'external'])
#         md_content = compose.Md(report=report,
#                                 header_file=header_file,
#                                 header='NYGC v7 runtime metrics',
#                                 section_header='Overview',
#                                 center_content=all_section_content,
#                                 figures=[])

        md_content = compose.Md(report=report,
                                header_file=header_file,
                                header='NYGC v7 runtime metrics',
                                section_header='Overview',
                                center_content=all_section_content,
                                figures=[(results.exit_status_table.script,
                                          results.exit_status_table.div),
                                          (results.disk_type_plot.script,
                                           results.disk_type_plot.div),
                                          (results.summary_table.script,
                                           results.summary_table.div)])
        #  =======================
        #  Task
        #  =======================
        section_content = ''.join([''])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        md_content.update_doc(section_header='Task metrics',
                              center_content=all_section_content,
                              figures=[(results.fig.script,
                                        results.fig.div)])
        #  =======================
        #  Subworkflow
        #  =======================
        section_content = ''.join([''])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        md_content.update_doc(section_header='Subworkflow metrics',
                              center_content=all_section_content,
                              figures=[(results.fig2.script,
                                        results.fig2.div)])
        #  =======================
        #  Workflow
        #  =======================
        section_content = ''.join([''])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        md_content.update_doc(section_header='Full pipeline metrics',
                              center_content=all_section_content,
                              figures=[(results.fig3.script,
                                        results.fig3.div)])

        
        
    
def main():
    output_info_file = sys.argv[1]
    metrics_file = sys.argv[2]
    plot_file = sys.argv[3]
    results = PlotRuntime(output_info_file=output_info_file,
                          metrics_file=metrics_file,
                          plot_file=plot_file)
    make_files(results, appendix=False)
    
    
if __name__ == "__main__":
    main()
         
        
