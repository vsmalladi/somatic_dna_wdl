import sys
import os
import json
import re
import logging as log
import pandas as pd
import numpy as np
import Colorer
import argparse
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
import plotly.graph_objects as go


pd.set_option('display.max_columns', 500)
log.basicConfig(format='%(levelname)s:  %(message)s', level=log.INFO)


class PlotRuntime():
    '''Plot results from outputMetrics.csv'''
    def __init__(self,
                 project_id,
                 manifest=False,
                 output_info_file=False,
                 metrics_file=False,
                 non_retry_metrics_file=False,
                 plot_file=False):
        self.set_colors()
        self.project_id = project_id
        if not plot_file:
            if manifest:
                basename = project_id.replace(' ', '_') + '_outputMetrics'
                self.out_md = basename + '.md'
                self.header_file = basename + '.header.txt'
                self.plot_file = basename + '.html'
            else:
                self.out_md = metrics_file.replace('.csv', '.md')
                self.header_file = metrics_file.replace('.csv', '.header.txt')
                self.plot_file = metrics_file.replace('.csv', '.html')
        else:
            self.out_md = plot_file.replace('.html', '.md')
            self.header_file = plot_file.replace('.html', '.header.txt')
            self.plot_file = plot_file
        if manifest:
            manifest_data = pd.read_csv(manifest)
            output_infos = []
            metadatas = []
            non_retry_metadatas = []
            for i, row in manifest_data.iterrows():
                output_info = self.load_input(row.output_info)
                metadata = self.load_metadata(row.metrics)
                non_retry_metadata = self.load_metadata(row.non_retried_metrics)
                non_retry_metadatas.append(non_retry_metadata)
                output_infos.append(output_info)
                metadatas.append(metadata)
            self.output_info = output_infos
            self.metadata = pd.concat(metadatas, ignore_index=True)
            self.non_retry_metadata = pd.concat(non_retry_metadatas, ignore_index=True)
        else:
            self.output_info = self.load_input(output_info_file)
            self.metadata = self.load_metadata(metrics_file)
            if os.path.isfile(non_retry_metrics_file):
            	self.non_retry_metadata = self.load_metadata(non_retry_metrics_file)
            else:
              	self.non_retry_metadata = pd.DataFrame({})
        # intro figs
        self.summary_table = self.get_table()
        self.gather_summary()
        if 'cpu_platform' in self.metadata.columns:
          self.plot_summary()
          self.plot_preempts()
          if 'wait_time_m' in self.metadata:
              self.plot_wait_time()
          if 'avg_total_cost' in self.metadata.columns:
              self.plot_cost()
          else:
              self.cost_fig = False
          # main figs
          self.fig = self.plot_by_steps(data_steps=self.metadata,
                                        non_retry_data_steps=self.non_retry_metadata,
                                        button=True,
                                        x="workflow_name",
                                        ys=["sample_task_core_h",
                                            "max_mem_used_gb",
                                            "sample_task_vm_runtime_h"],
                                        color="task_call_name",
                                        hover_data=['id'],
                                        color_discrete_sequence=self.colors_set2,
                                        title='<b>Tasks: resource usage summary plot</b>',
                                        labels={'sample_task_vm_runtime_h' : 'Task runtime(h)',
                                                'max_mem_used_gb' : 'Mem (G)',
                                                'sample_task_core_h' : 'Task core hours',
                                                'workflow_name' : 'Sub-workflow',
                                                'task_call_name' : 'Task'},
                                        x_title="Pipeline subworkflows")
          self.fig2 = self.plot_by_steps(data_steps=self.metadata,
                                         non_retry_data_steps=self.non_retry_metadata,
                                         button=False,
                                         x="workflow_name",
                                         ys=["sample_subworkflow_core_h",
                                            "subworkflow_max_mem_g",
                                            "sample_subworkflow_run_time_h"],
                                         color="workflow_name",
                                         hover_data=['id'],
                                         color_discrete_sequence=self.colors_set1,
                                         labels={'sample_subworkflow_run_time_h' : 'sub-workflow runtime(h)',
                                                'subworkflow_max_mem_g' : 'Mem (G)',
                                                'sample_subworkflow_core_h' : 'Sub-workflow core hours',
                                                'workflow_name' : 'Sub-workflow',
                                                'task_call_name' : 'Task'},
                                         title='<b>Subworkflows: resource usage summary plot</b>',
                                         x_title="Pipeline subworkflows")
          self.fig3 = self.plot_by_steps(data_steps=self.metadata,
                                         non_retry_data_steps=self.non_retry_metadata,
                                         button=False,
                                         x="id",
                                         ys=["sample_workflow_core_h",
                                            "workflow_max_mem_g",
                                            "sample_workflow_run_time_h"],
                                         color="id",
                                         hover_data=['id'],
                                         color_discrete_sequence=self.colors,
                                         labels={'sample_workflow_run_time_h' : 'workflow runtime(h)',
                                                'workflow_max_mem_g' : 'Mem (G)',
                                                'sample_workflow_core_h' : 'Workflow core hours',
                                                'workflow_name' : 'Sub-workflow',
                                                'task_call_name' : 'Task'},
                                         title='<b>Workflows: resource usage summary plot</b>',
                                         x_title="Full pipeline")
          self.fig4, self.fig5, self.fig6, self.cummulative_cost_fig = self.custom_by_steps(data_steps=self.metadata,
                                                                non_retry_data_steps=self.non_retry_metadata,
                                                                button=True,
                                                                x="workflow_name",
                                                                top_ys=["mem_total_gb", "disk_total_gb", ""],
                                                                central_ys=["max_mem_used_gb", "max_disk_used_gb", "sample_task_vm_runtime_h"],
                                                                color="task_call_name",
                                                                hover_data=['id'],
                                                                color_discrete_sequence=self.colors,
                                                                title='<b>Tasks: resources used vs available</b>',
                                                                labels={'max_disk_used_gb' : 'Disk used (G)',
                                                                        'disk_total_gb' : 'Available disk (G)',
                                                                        'mem_total_gb' : 'Available mem (G)',
                                                                        'sample_task_vm_runtime_h' : 'Task runtime(h)',
                                                                        'max_mem_used_gb' : 'Mem (G)',
                                                                        'sample_task_core_h' : 'Task core hours',
                                                                        'workflow_name' : 'Sub-workflow',
                                                                        'task_call_name' : 'Task'},
                                                                x_title="Pipeline subworkflow")
        else:
          log.warning('Monitoring_image not used not enough runtime metrics to create plot')
          sys.exit(0)
    
        
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
        self.metadata = self.metadata[~(self.metadata.backend_status.isin(['Running', 'CacheHit']))].copy()
        # exit status of non-zero exits
        exit_status = self.metadata[['task_call_name',
                                     'backend_status']].value_counts().to_frame().reset_index()
        self.exit_status = exit_status[exit_status.backend_status != 'Success']
        self.exit_status.columns = ['Task', 'Exit status', 'Count']
        # preemptible percentage
        testable = self.metadata.dropna(subset=['preemptible'])
        self.preemptible_percent = (testable[testable.preemptible].shape[0] / float(testable.shape[0])) * 100
        self.preemptible_percent_runtime = (testable[testable.preemptible].vm_runtime_m.sum() / float(testable.vm_runtime_m.sum())) * 100
        # disk type
        self.metadata['disk_type']= self.metadata.apply(lambda row: row.disk_types.replace('[\'', '').replace('\']', '')
                                                        if str(row.disk_types) != 'nan' else '', axis=1)
        # create vm summary
        if 'machine_type' in self.metadata.columns and 'cpu_platform' in self.metadata.columns:
            self.metadata['vm_summary'] = self.metadata.cpu_platform + ' ' + self.metadata.machine_type + ' ' + self.metadata.disk_type
        elif 'cpu_platform' in self.metadata.columns:
          	self.metadata['vm_summary'] = self.metadata.cpu_platform + ' ' + self.metadata.disk_type
        else:
            self.metadata['vm_summary'] = self.metadata.disk_type
        self.disk_type = self.metadata.groupby(['task_call_name', 'id', 'instance_name', 'disk_type', 'vm_summary']).sample_task_vm_runtime_h.sum().reset_index()
        self.disk_type.columns = ['Task', 'Id', 'Instance name', 'Disk type', 'VM type', 'Wall clock (h)']
        
    def add_button(self, fig):
        fig.update_layout(updatemenus=[{'type' : 'buttons',
                                        'direction' : 'up',
                                        'buttons' : [{'label' : 'All instances',
                                                    'method' : 'update',
                                                    'args' : [{'visible' : [True] * self.all_levels + [False] * self.non_retry_levels},
                                                              {'title' : 'All instances'}]
                                                    },
                                                    {'label' : 'Non-Preempted instances',
                                                    'method' : 'update',
                                                    'args' : [{'visible' : [False] * self.all_levels + [True] * self.non_retry_levels},
                                                              {'title' : 'Non-Preempted instances'}]
                                                    }]
                                       }])
        return fig
    
    def plot_preempts(self):
        preempts = self.metadata[['task_call_name', 'execution_status', 'sample_task_vm_runtime_h']][self.metadata['execution_status'] != 'Done'].copy()
        if not preempts.empty:
            preempts.columns = ['Task', 'Exit status', 'Wall clock (h)']
            fig = px.box(preempts, x='Task', y='Wall clock (h)', points="all",
                              color_discrete_sequence=self.colors_set2,
                              color='Exit status')
            fig.update_layout(xaxis_type='category',
                              title_text='Runtime of preempted or failed jobs')
            fig.update_xaxes(title_text='')
            fig.update_yaxes(title_text="Wall clock (h)")
            self.preempt_fig = ploter.Fig(fig)
        else:
            self.preempt_fig = False
    
    def plot_wait_time(self):
        waits = self.metadata[['task_call_name', 'vm_summary', 'sample_task_vm_runtime_h', 'wait_time_m']].copy()
        if not waits.empty:
            waits['vm_summary'].fillna('Unknown', inplace=True)
            waits = waits.dropna(subset=['sample_task_vm_runtime_h', 'wait_time_m'])
            waits.columns = ['Task', 'VM type', 'Wall clock (h)', 'Wait time (m)']
            fig = px.box(waits, x='Task', y='Wait time (m)', points="all",
                              color_discrete_sequence=self.colors_set1,
                              color='VM type')
            fig.update_layout(xaxis_type='category',
                              title_text='Wait time before jobs begin (m)')
            fig.update_xaxes(title_text='')
            fig.update_yaxes(title_text="Wait time (m)")
            self.wait_fig = ploter.Fig(fig)
        else:
            self.wait_fig = False
    
    def plot_cost(self):
        costs = self.metadata[['task_call_name', 'vm_summary', 'sample_task_vm_runtime_h',
                                'avg_total_cost', 'runtime_scaled_total_cost']].copy()
        if not costs.empty:
            costs['vm_summary'].fillna('Unknown', inplace=True)
            costs = costs.dropna(subset=['sample_task_vm_runtime_h', 'runtime_scaled_total_cost'])
#            costs = costs.dropna(subset=['runtime_scaled_total_cost']).copy()
            costs.columns = ['Task', 'VM type', 'Wall clock (h)', 'Average cost', 'Runtime scaled cost']
            fig1 = px.box(costs, x='Task', y='Runtime scaled cost', points="all",
                              color_discrete_sequence=self.colors_set1,
                              color='VM type')
            fig1.update_layout(xaxis_type='category',
                              title_text='Cost')
            fig1.update_xaxes(title_text='')
            fig1.update_yaxes(title_text="Runtime scaled cost")
            # create scatter plot
            fig2 = px.scatter(costs, x='Wall clock (h)', y='Runtime scaled cost',
                              color_discrete_sequence=self.colors_set2,
                              color='Task')
            fig2.update_layout(title_text='Cost by runtime')
            fig2.update_xaxes(title_text='Wall clock (h)')
            fig2.update_yaxes(title_text="Runtime scaled cost")
            fig = make_subplots(rows=1,
                        cols=2,
                        subplot_titles=("", ""))
            # All instances
            for trace in fig1.data:
                trace['showlegend'] = True
                trace['pointpos'] = 0
                fig.add_trace(
                              trace,
                              row=1, col=1
                             )
            for trace in fig2.data:
                trace['showlegend'] = True
                fig.add_trace(
                              trace,
                              row=1, col=2
                           )
            # add labels
            fig.update_layout(title_text="Cost details")
            fig.update_xaxes(title_text='')
            fig.update_yaxes(title_text="Runtime-scaled cost",
                             range=[0, costs['Runtime scaled cost'].max() + 2],
                             row=1, col=1)
            fig.update_yaxes(title_text="",
                             range=[0, costs['Runtime scaled cost'].max() + 2],
                             row=1, col=2)
            fig.update_xaxes(title_text="",
    #                         categoryarray=category_orders[ys[0]],
    #                         categoryorder='array',
                             row=1, col=1)
            fig.update_xaxes(title_text="Wall clock (h)",
                             row=1, col=2)
    #        fig.update_xaxes(tickfont={'size' : 5})
            self.cost_fig = ploter.Fig(fig)
        else:
            self.cost_fig = False
        
    def plot_summary(self):
        fig = px.box(self.disk_type, x='VM type', y='Wall clock (h)', points="all",
                          color_discrete_sequence=self.colors_set2,
                          color='Task',
                          height=300)
        for trace in fig.data:
            trace['showlegend'] = True
            trace['pointpos'] = 0
        fig.update_layout(xaxis_type='category',
                          title_text='<b>Runtime by VM type</b>')
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
        lines = ['Below is a table of all the tasks with non-zero exit status. ']
        lines += ['Overall, ' + str(int(self.preemptible_percent)) + '% of instances were preemptible. ',
                  '' + str(int(self.preemptible_percent_runtime)) + '% of total run time was on preemptible instances.']
        return lines

    def get_table(self):
        '''Make workflow resource usage summary table'''
        table_data = self.metadata[['id', 'sample_workflow_run_time_h']].copy()
        table_data = table_data.dropna(subset=['sample_workflow_run_time_h'])
        table_data = table_data.sort_values(['sample_workflow_run_time_h']).copy()
        table_data.columns = ['Id', 'Wallclock (h)']
        summary_table = ploter.make_html_table(table_data.drop_duplicates(),
                                               name='runtime_overview',
                                               col_widths=['80%','20%'],
                                               searching=True,
                                               scroll=True)
        return summary_table
    
    def custom_by_steps(self,
                          data_steps,
                          non_retry_data_steps,
                          button=False,
                          x="workflow_name",
                          top_ys=["mem_total_gb", "disk_total_gb", ""],
                          central_ys=["max_mem_used_gb", "max_disk_used_gb", "sample_task_vm_runtime_h"],
                          color="task_call_name",
                          hover_data=['id'],
                          color_discrete_sequence=False,
                          title='',
                          labels={},
                          x_title="Pipeline subworkflow"):
        ''' Plot Task disk used vs total dis, Mem vs max mem and stacked barplot of total time taken'''
        category_orders = {}
        color_cols = [col for col in [x, color] + top_ys + central_ys if not col == '']
        data_steps = data_steps.drop_duplicates(subset=color_cols).copy()
        data_steps['task per id'] = data_steps.apply(lambda row: ' '.join([row[x], row['id']]), axis=1)
        for y in central_ys:
            if y == 'sample_task_vm_runtime_h':
                grouped = data_steps.groupby(['task per id'])
                category_orders[y] = grouped[y].agg(['sum']).reset_index().sort_values('sum')['task per id'].unique().tolist()
            else:
                grouped = data_steps.groupby(x)
                category_orders[y] = grouped[y].agg(['max']).reset_index().sort_values('max')[x].unique().tolist()
        # plot total vs max mem
        max_color = '#C0FF02'
        max_colors = [max_color] * data_steps.shape[0]
        fig1 = px.box(data_steps, x=x, y=central_ys[0], points="all",
                      hover_data=['id', top_ys[0], central_ys[0]],
                      color_discrete_sequence=color_discrete_sequence,
                      labels=labels,
                      boxmode='overlay',
                      color=color, category_orders=category_orders)
        for trace in fig1.data:
            trace['pointpos'] = 0
        fig1.add_trace(go.Scatter(
            x=data_steps[x], y=data_steps[top_ys[0]],
            text=data_steps['task_call_name'],
            mode='markers',
            name='Mem requested',
            marker_color=max_color
            ))
        # plot total vs max disk used
        fig2 = px.box(data_steps, x=x, y=central_ys[1], points="all",
                      hover_data=['id', top_ys[1], central_ys[1]],
                      labels=labels,
                      boxmode='overlay',
                      color_discrete_sequence=color_discrete_sequence,
                      color=color, category_orders=category_orders)
        for trace in fig2.data:
            trace['pointpos'] = 0
        fig2.add_trace(go.Scatter(
            x=data_steps[x], y=data_steps[top_ys[1]],
            text=data_steps['task_call_name'],
            mode='markers',
            name='Disk requested',
            marker_color=max_color
            ))
        # plot stacked bar of task persample summed runtime [x, 'id']
        data_steps['Result set'] = 'All runs'
        non_preempted_data_steps = data_steps[data_steps.backend_status == 'Success'].copy()
        non_preempted_data_steps['Result set'] = 'Non-preemptedruns'
        plot_data = pd.concat([data_steps, non_preempted_data_steps], ignore_index=True)
        fig3 = px.bar(plot_data, x='task per id', y=central_ys[2],
                      hover_data=['id', 'task per id', central_ys[2], 'preemptible'],
                      labels=labels,
                      facet_row='Result set',
                      title='<b>Cumulative runtime per ID</b>',
                      color_discrete_sequence=color_discrete_sequence,
                      color=color, category_orders=category_orders)
        fig3.update_traces(marker=dict(line=dict(width=0)))
        fig3.update_xaxes(tickangle=45)
        fig3.update_xaxes(matches='x')
        if x.lower() in ['task per id', 'id']:
            fig.update_xaxes(tickfont={'size' : 5})
        if 'avg_total_cost' in data_steps.columns:
            # all tasks
            fig4 = px.bar(data_steps, x='task per id', y='runtime_scaled_total_cost',
                          hover_data=['id', 'task per id', 'runtime_scaled_total_cost', 'preemptible'],
                          labels=labels,
                          facet_row='Result set',
                          title='<b>Cumulative cost per ID and subworkflow</b>',
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
            fig4.update_traces(marker=dict(line=dict(width=0)))
            fig4.update_xaxes(tickangle=45)
            fig4.update_xaxes(matches='x')
            if x.lower() in ['task per id', 'id']:
                fig.update_xaxes(tickfont={'size' : 5})
            # succeeded tasks only
            fig5 = px.bar(non_preempted_data_steps, x='task per id', y='runtime_scaled_total_cost',
                          hover_data=['id', 'task per id', 'runtime_scaled_total_cost', 'preemptible'],
                          labels=labels,
                          facet_row='Result set',
                          title='<b>Cumulative cost per ID</b>',
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
            fig5.update_traces(marker=dict(line=dict(width=0)))
            fig5.update_xaxes(tickangle=45)
            fig5.update_xaxes(matches='x')
            if x.lower() in ['task per id', 'id']:
                fig.update_xaxes(tickfont={'size' : 5})
            fig = make_subplots(rows=1,
                    cols=1,
                    subplot_titles=(""))
            # All instances
            for trace in fig4.data:
                trace['showlegend'] = True
                fig.add_trace(
                              trace,
                              row=1, col=1
                             )
            for trace in fig5.data:
                trace['showlegend'] = True
                trace['visible'] = False
                fig.add_trace(
                              trace,
                              row=1, col=1
                           )
#            if button:
            self.all_levels = len(fig4.data)
            self.non_retry_levels = len(fig5.data)
            fig = self.add_button(fig)
            fig.update_layout(barmode='stack')
            fig.update_layout(title_text="Cumulative cost per ID and subworkflow")
            fig.update_xaxes(title_text='')
            data_steps['summed_cost'] = data_steps.groupby(['task per id'])['runtime_scaled_total_cost'].transform('sum')
            fig.update_yaxes(title_text="Runtime-scaled cost",
                             range=[0, data_steps['summed_cost'].max() + 2],
                             row=1, col=1)
            cummulative_cost_fig = ploter.Fig(fig)
        else:
            cummulative_cost_fig = False
        return ploter.Fig(fig1), ploter.Fig(fig2), ploter.Fig(fig3), cummulative_cost_fig
    
    def plot_by_steps(self,
                      data_steps,
                      non_retry_data_steps,
                      button=False,
                      x="workflow_name",
                      ys=["sample_subworkflow_core_h", "mem_total_gb", "sample_task_vm_runtime_h"],
                      color="task_call_name",
                      hover_data=['id'],
                      color_discrete_sequence=False,
                      title='',
                      labels={},
                      x_title="Pipeline subworkflow"):
        ''' Plot Task CPU, Mem and Wallclock'''
        category_orders = {}
        data_steps = data_steps.drop_duplicates(subset=[x, color] + ys).copy()
        for y in ys:
            grouped = data_steps.groupby(x)
            category_orders[y] = grouped[y].agg(['max']).reset_index().sort_values('max')[x].unique().tolist()
        # All instances
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
        self.all_levels = len(fig1.data) + len(fig2.data) + len(fig3.data)
        if non_retry_data_steps.empty:
          button = False
        if button:
            # Non-RetryableFailure instances
            fig4 = px.box(non_retry_data_steps, x=x, y=ys[0], points="all",
                          hover_data=hover_data,
                          color_discrete_sequence=color_discrete_sequence,
                          labels=labels,
                          color=color, category_orders=category_orders)
            fig5 = px.box(non_retry_data_steps, x=x, y=ys[1], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
            fig6 = px.box(non_retry_data_steps, x=x, y=ys[2], points="all",
                          hover_data=hover_data,
                          labels=labels,
                          color_discrete_sequence=color_discrete_sequence,
                          color=color, category_orders=category_orders)
            self.non_retry_levels = len(fig4.data) + len(fig5.data) + len(fig6.data)
        fig = make_subplots(rows=1,
                            cols=3,
                            subplot_titles=("", " ", " "))
        # All instances
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
        # Non-Preempted instances
        if button:
            for trace in fig4.data:
                trace['showlegend'] = False
                trace['visible'] = False
                trace['pointpos'] = 0
                fig.add_trace(
                              trace,
                              row=1, col=1
                             )
            for trace in fig5.data:
                trace['showlegend'] = False
                trace['visible'] = False
                trace['pointpos'] = 0
                fig.add_trace(
                              trace,
                              row=1, col=2
                           )
            for trace in fig6.data:
                trace['showlegend'] = True
                trace['visible'] = False
                trace['pointpos'] = 0
                fig.add_trace(
                              trace,
                              row=1, col=3
                           )
        fig.update_layout(xaxis_type='category',
                          title_text=title)
        fig.update_xaxes(title_text='')
        fig.update_yaxes(title_text="CPU Time (h)",
                         range=[0, self.metadata[ys[0]].max() + 2],
                         row=1, col=1)
        fig.update_yaxes(title_text="Mem (G)",
                         range=[0, self.metadata[ys[1]].max() + 2],
                         row=1, col=2)
        fig.update_yaxes(title_text="Wall clock (h)",
                         range=[0, self.metadata[ys[2]].max() + 2],
                         row=1, col=3)
        fig.update_xaxes(title_text="",
                         categoryarray=category_orders[ys[0]],
                         categoryorder='array',
                         row=1, col=1)
        fig.update_xaxes(title_text="",
                         categoryarray=category_orders[ys[1]],
                         categoryorder='array',
                         row=1, col=2)
        fig.update_xaxes(title_text="",
                         categoryarray=category_orders[ys[2]],
                         categoryorder='array',
                         row=1, col=3)
        if x.lower() == 'id':
            fig.update_xaxes(tickfont={'size' : 5})
        if button:
            fig = self.add_button(fig)
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
#         basic_info = ''.join(['\n\nRuntime metrics: v7 pipeline',
#                               '\n\nRuntime for different task, subworkflows and full workflows of the v7 pipeline are shown below.\n\n',
#                               '\n\nTask runtime and core hours can have preempted workflow runtime removed but this cannot be subtracted from the full workflow runtime.\n\n'
#                               ])
        basic_info = ''.join(['\n\n<b>Project :</b> ', results.project_id, '\n\n'])
        logo = '<img src="NY-Genome-Logo.jpg" class="center">' + ('\n' * 4)
        section_content = ''.join([logo, basic_info, '\n\n'])
        all_section_content = contents.add_line(content=section_content,
                                               libraries=['WGS', 'Exome'],
                                               types=['tumor_only', 'paired'],
                                               audiences=['internal', 'external'])
        md_content = compose.Md(report=report,
                                header_file=header_file,
                                header= results.project_id + ' NYGC v7',
                                section_header='Project overview',
                                center_content=all_section_content,
                                figures=[])
        #  =======================
        #  Runtime
        #  =======================
        section_content = ''.join([''])
#         section_content = ''.join(['End-to-end runtime per sample'])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        md_content.update_doc(section_header='Runtime overview',
                              center_content=all_section_content,
                              figures=[(results.summary_table.script,
                                           results.summary_table.div)])
        #  =======================
        #  Non-zero exit status
        #  =======================
        if results.preempt_fig:
            section_content = ''.join(results.write_summary())
    #         section_content = ''.join(['All tasks with non-zero exit status.'])
            lines = [contents.add_line(content=section_content,
                                        libraries=['WGS', 'Exome'],
                                        types=['tumor_only', 'paired'],
                                        audiences=['internal', 'external'])]
            all_section_content = md_content.join_ignore_none(lines)
            md_content.update_doc(section_header='Non-zero exits',
                                  center_content=all_section_content,
                                  figures=[(results.exit_status_table.script,
                                              results.exit_status_table.div),
                                            (results.preempt_fig.script,
                                             results.preempt_fig.div)])
        #  =======================
        #  Wait time (m)
        #  =======================
        if results.wait_fig:
            section_content = ''.join([''])
    #         section_content = ''.join(['All tasks with non-zero exit status.'])
            lines = [contents.add_line(content=section_content,
                                        libraries=['WGS', 'Exome'],
                                        types=['tumor_only', 'paired'],
                                        audiences=['internal', 'external'])]
            all_section_content = md_content.join_ignore_none(lines)
            md_content.update_doc(section_header='Job wait time',
                                  center_content=all_section_content,
                                  figures=[(results.wait_fig.script,
                                            results.wait_fig.div)])
        #  =======================
        #  Cost
        #  =======================
        if results.cost_fig:
            section_content = ''.join([''])
    #         section_content = ''.join(['All tasks with non-zero exit status.'])
            lines = [contents.add_line(content=section_content,
                                        libraries=['WGS', 'Exome'],
                                        types=['tumor_only', 'paired'],
                                        audiences=['internal', 'external'])]
            all_section_content = md_content.join_ignore_none(lines)
            md_content.update_doc(section_header='Cost overview',
                                  center_content=all_section_content,
                                  figures=[(results.cost_fig.script,
                                            results.cost_fig.div),
                                            (results.cummulative_cost_fig.script,
                                            results.cummulative_cost_fig.div)])
        #  =======================
        #  Task
        #  =======================
        section_content = ''.join([''])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        figures=[(results.fig.script,
                results.fig.div),
                (results.fig4.script,
                results.fig4.div),
                (results.fig5.script,
                results.fig5.div),
                (results.fig6.script,
                results.fig6.div)]
        md_content.update_doc(section_header='Task metrics',
                              center_content=all_section_content,
                              figures=figures)
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
        #  =======================
        #  Disk types
        #  =======================
        section_content = ''.join(['Disk type by runtime for each task in the run.'])
        lines = [contents.add_line(content=section_content,
                                    libraries=['WGS', 'Exome'],
                                    types=['tumor_only', 'paired'],
                                    audiences=['internal', 'external'])]
        all_section_content = md_content.join_ignore_none(lines)
        md_content.update_doc(section_header='Disk type',
                              center_content=all_section_content,
                              figures=[(results.disk_type_plot.script,
                                           results.disk_type_plot.div)])

        
def get_args():
    '''Parse input flags
    '''
    parser = argparse.ArgumentParser()
    parser.add_argument('--output-info',
                        help='JSON file with outputInfo. Includes workflow uuid, '
                        'options JSON dictionary submitted to cromwell, '
                        'project pairing, sample, '
                        'genome build, library and interval list information. '
                        'Also includes output files and '
                        'the sub-workflow uuids',
                        required=False
                        )
    parser.add_argument('--plot',
                        help='Output path for html file.',
                        required=True
                        )
    parser.add_argument('--metrics',
                        help='CSV file with outputMetrics. '
                        'Includes sample id mem, wallclock time, core hours '
                        'for tasks with any exit status (including Failed)',
                        required=False
                        )
    parser.add_argument('--non-retry-metrics',
                        help='CSV file with outputMetrics.non_retried. '
                        'Includes sample id mem, wallclock time, core hours '
                        'for tasks with an exit status of zero',
                        required=False
                        )
    parser.add_argument('--multi',
                        help='Records for multiple runs are stored in the JSON file with outputInfo.',
                        required=False,
                        action='store_true'
                        )
    parser.add_argument('--name',
                        default='nygc_pipeline',
                        help='Use with the --multi flag to name the output.',
                        required=False
                        )
    parser.add_argument('--manifest',
                        help='Use with the --multi flag to list outputInfo JSON files '
                        'outputMetrics.non_retried files and outputMetrics files.',
                        required=False
                        )
    args_namespace = parser.parse_args()
    return args_namespace.__dict__
    
def main():
    args = get_args()
    if args['multi']:
        results = PlotRuntime(project_id=args['name'],
                              manifest=args['manifest'],
                              plot_file=args['plot'])
    else:
        results = PlotRuntime(project_id=args['name'],
                              output_info_file=args['output_info'],
                              metrics_file=args['metrics'],
                              non_retry_metrics_file=args['non_retry_metrics'],
                              plot_file=args['plot'])
    make_files(results, appendix=False)
    
    
if __name__ == "__main__":
    main()
         
        
