import plotly.offline as po
import pandas as pd
import re
from bs4 import BeautifulSoup


class Fig():
    ''' Create the script and div string for a single figure'''

    def __init__(self, fig,
                 scale=False,
                 automargin=True,
                 remove_buttons=['sendDataToCloud',
                                 'toggleSpikelines',
                                 'autoScale2d',
                                 'select2d',
                                 'lasso2d']):
        self.scale = scale
        self.fig = self.style_plot(fig)
        self.fig = self.fix_legend(self.fig)
        self.fig.update_yaxes(automargin=automargin)
        self.fig.update_yaxes(showgrid=False, zeroline=False)
        self.fig.update_xaxes(showgrid=False, zeroline=False)
        for i, annotation in enumerate(self.fig.layout['annotations']):
            self.fig.layout['annotations'][i]['font'] = {'color' : '#0080ff',
                                                            'size' : 8
                                                        }
        if self.scale:
            self.fig = self.scale_sizes(self.fig)
        self.remove_buttons = remove_buttons
        self.script, self.div = self.create_html()

    def style_plot(self, fig):
        '''
        Set standard layout.
        '''
        fig.update_layout(template='plotly_white', margin=dict(l=5,r=5,t=25,b=15), font=dict(size=8, color='#0080ff'))
        return fig
    
    def scale_sizes(self, fig):
        for i, data in enumerate(fig.data):
            new_sizes = [(size + 400) for size in data['marker']['size']]
            fig.data[i]['marker']['size'] = new_sizes
        return fig

    def fix_legend(self, fig):
        for i, data in enumerate(fig.data):
            try:
                fig.data[i].legendgroup = fig.data[i].legendgroup.split('=')[-1] # legend name and clear hovername
                fig.data[i].name = fig.data[i].name.split('=')[-1]
            except AttributeError:
                pass
        return fig
    
    def clean_tag(self, tag):
        '''clean up tag'''
        return re.sub(r' \s+', ' ', str(tag))

    def create_html(self):
        '''
         Get div for plotly fig
         ['sendDataToCloud', 'toggleSpikelines', 'hoverCompareCartesian',
         'zoom2d', 'zoomIn2d', 'zoomOut2d',
         'autoScale2d', 'resetScale2d', 'select2d',
         'lasso2d']
         '''
        scripts = []
        divs = []
        plotly_div = po.plot(self.fig,
                             include_plotlyjs=False,
                             output_type='div',
                             config={'modeBarButtonsToRemove':
                                     self.remove_buttons},
                             auto_open=False)
        soup = BeautifulSoup(plotly_div, "html.parser")
        for tag in soup.find_all():
            if tag.name == 'script':
                scripts.append(tag)
            elif tag.name == 'div':
                divs.append(tag)
        script = '\n' + '\n'.join([self.clean_tag(tag) for tag in scripts]) + '\n'
        div = '\n' + '\n'.join([self.clean_tag(tag) for tag in divs]) + '\n'
        return script, div
    
    
def make_html_table(table_data, name, col_widths, searching=False, scroll=False):
    
    pd.set_option('display.max_colwidth', None)
    if len(table_data.columns) != len(col_widths):
        log.error('Number of columns in table and col_widths does not match')
        return
    class html_table:
        width_string = ''
        for i, width in enumerate(col_widths):
            width_string = width_string + '{"width": "'+ width +'"},'
        if searching:
            search_bool = 'true'
        else:
            search_bool = 'false'

        if scroll:
            scroll_string = '"scrollY": "320px", "scrollCollapse": true,'
        else:
            scroll_string = ''

        script = '''<script>
            $(document).ready(function() {
                $("#''' + name + '''").DataTable({
                    dom: 'Bfrtip',
                    "columns": [''' + width_string + '''],
                    "paging": false,
                    "searching":''' + search_bool + ''',
                    "info": false,
                    "autoWidth": false,
                    buttons: [{extend:'copy', text: 'Copy to clipboard'}],
                    "order": [],
                    ''' + scroll_string + '''
                });
                $($.fn.dataTable.tables( true ) ).DataTable().columns.adjust().draw();
                $(window).resize(function(){
                    $($.fn.dataTable.tables( true ) ).DataTable().columns.adjust().draw();
                });

            });
        </script>'''
        div = table_data.to_html(table_id=name, index=False, classes=['display '])
    return html_table
