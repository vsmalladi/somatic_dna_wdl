import os
import re

# docs
import Colorer
import logging as log
log.basicConfig(level=log.DEBUG)

class Content(object):
    '''
        Store lines and metadata about which lines should be included
        based on project metadata.
    '''
    def __init__(self,
                 tumor_only,
                 internal,
                 library):
        # set doc top
        self.library = library
        if tumor_only:
            self.type = 'tumor_only' # tumor-only or paired
        else:
            self.type = 'paired'
        if internal:
            self.audience = 'internal' # tumor-only or paired
        else:
            self.audience = 'external'
    
    def add_line(self,
                 content,
                 libraries=['WGS', 'Exome'],
                 types=['tumor_only', 'paired'],
                 audiences=['internal', 'external']):
        self.libraries = libraries
        self.types = types
        self.audiences = audiences
        self.content = content
        return self.check_lines()
    
    def check_lines(self):
        '''
        assemble text from
        '''
        if self.library in self.libraries \
                and self.type in self.types \
                and self.audience in self.audiences:
            return self.content


class Md(object):
    '''
        Write MD for pandoc and include HTML plotly plots.
    '''

    def __init__(self, report,
                 header_file,
                 header='NYGC Annotated Somatic Project Report',
                 default='show',
                 section_header='Overview',
                 center_content='''Section description''',
                 figures = None,
                 format='html'):
        # set doc top
        self.header_file = header_file
        self.default = default
        self.header = ['%' + header]
        # '<link rel="stylesheet" href="https://maxcdn.bootstrapcdn.com/bootstrap/3.3.1/css/bootstrap.min.css">',
        new = ['<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>']
        self.header_includes = ['<link rel="icon" type="image/png" href="favicon.png"/>',
                                '<script type="text/javascript" src="https://code.jquery.com/jquery-3.3.1.js"></script>',
                                '<script type="text/javascript" src="https://cdn.datatables.net/v/bs/dt-1.10.18/b-1.5.6/b-html5-1.5.6/datatables.min.js"></script>',
                                '''<style type="text/css">

                                    table {
                                        font-size:11px;
                                    }
                                    thead {
                                        color: white;
                                        background: dodgerblue;
                                    }
                                    #hla_table {
                                        width: 40% !important;
                                    }
                                    #mantis_table {
                                        width: 40% !important;
                                    }
                                    #mut_sigs_table {
                                        width: 40% !important;
                                    }
                                    .form-control {
                                        background-color: white;
                                    }
                                    /*
                                     * This combined file was created by the DataTables downloader builder:
                                     *   https://datatables.net/download
                                     *
                                     * To rebuild or modify this file with the latest versions of the included
                                     * software please visit:
                                     *   https://datatables.net/download/#bs/dt-1.10.18/b-1.5.6/b-html5-1.5.6
                                     *
                                     * Included libraries:
                                     *   DataTables 1.10.18, Buttons 1.5.6, HTML5 export 1.5.6
                                     */

                                    table.dataTable {
                                      clear: both;
                                      margin-top: 6px !important;
                                      margin-bottom: 6px !important;
                                      max-width: none !important;
                                      border-collapse: separate !important;
                                    }
                                    table.dataTable td,
                                    table.dataTable th {
                                      -webkit-box-sizing: content-box;
                                      box-sizing: content-box;
                                    }
                                    table.dataTable td.dataTables_empty,
                                    table.dataTable th.dataTables_empty {
                                      text-align: center;
                                    }
                                    table.dataTable.nowrap th,
                                    table.dataTable.nowrap td {
                                      white-space: nowrap;
                                    }

                                    div.dataTables_wrapper div.dataTables_length label {
                                      font-weight: normal;
                                      text-align: left;
                                      white-space: nowrap;
                                    }
                                    div.dataTables_wrapper div.dataTables_length select {
                                      width: 75px;
                                      display: inline-block;
                                    }
                                    div.dataTables_wrapper div.dataTables_filter {
                                      text-align: right;
                                    }
                                    div.dataTables_wrapper div.dataTables_filter label {
                                      font-weight: normal;
                                      white-space: nowrap;
                                      text-align: left;
                                    }
                                    div.dataTables_wrapper div.dataTables_filter input {
                                      margin-left: 0.5em;
                                      display: inline-block;
                                      width: auto;
                                    }
                                    div.dataTables_wrapper div.dataTables_info {
                                      padding-top: 8px;
                                      white-space: nowrap;
                                    }
                                    div.dataTables_wrapper div.dataTables_paginate {
                                      margin: 0;
                                      white-space: nowrap;
                                      text-align: right;
                                    }
                                    div.dataTables_wrapper div.dataTables_paginate ul.pagination {
                                      margin: 2px 0;
                                      white-space: nowrap;
                                    }
                                    div.dataTables_wrapper div.dataTables_processing {
                                      position: absolute;
                                      top: 50%;
                                      left: 50%;
                                      width: 200px;
                                      margin-left: -100px;
                                      margin-top: -26px;
                                      text-align: center;
                                      padding: 1em 0;
                                    }

                                    table.dataTable thead > tr > th.sorting_asc, table.dataTable thead > tr > th.sorting_desc, table.dataTable thead > tr > th.sorting,
                                    table.dataTable thead > tr > td.sorting_asc,
                                    table.dataTable thead > tr > td.sorting_desc,
                                    table.dataTable thead > tr > td.sorting {
                                      padding-right: 30px;
                                    }
                                    table.dataTable thead > tr > th:active,
                                    table.dataTable thead > tr > td:active {
                                      outline: none;
                                    }
                                    table.dataTable thead .sorting,
                                    table.dataTable thead .sorting_asc,
                                    table.dataTable thead .sorting_desc,
                                    table.dataTable thead .sorting_asc_disabled,
                                    table.dataTable thead .sorting_desc_disabled {
                                      cursor: pointer;
                                      position: relative;
                                    }
                                    table.dataTable thead .sorting:after,
                                    table.dataTable thead .sorting_asc:after,
                                    table.dataTable thead .sorting_desc:after,
                                    table.dataTable thead .sorting_asc_disabled:after,
                                    table.dataTable thead .sorting_desc_disabled:after {
                                      position: absolute;
                                      bottom: 8px;
                                      right: 8px;
                                      display: block;
                                      font-family: 'Glyphicons Halflings';
                                      opacity: 0.5;
                                    }
                                    table.dataTable thead .sorting:after {
                                      opacity: 0.2;
                                      content: "\e150";
                                      /* sort */
                                    }
                                    table.dataTable thead .sorting_asc:after {
                                      content: "\e155";
                                      /* sort-by-attributes */
                                    }
                                    table.dataTable thead .sorting_desc:after {
                                      content: "\e156";
                                      /* sort-by-attributes-alt */
                                    }
                                    table.dataTable thead .sorting_asc_disabled:after,
                                    table.dataTable thead .sorting_desc_disabled:after {
                                      color: #eee;
                                    }

                                    div.dataTables_scrollHead table.dataTable {
                                      margin-bottom: 0 !important;
                                    }

                                    div.dataTables_scrollBody > table {
                                      border-top: none;
                                      margin-top: 0 !important;
                                      margin-bottom: 0 !important;
                                    }
                                    div.dataTables_scrollBody > table > thead .sorting:after,
                                    div.dataTables_scrollBody > table > thead .sorting_asc:after,
                                    div.dataTables_scrollBody > table > thead .sorting_desc:after {
                                      display: none;
                                    }
                                    div.dataTables_scrollBody > table > tbody > tr:first-child > th,
                                    div.dataTables_scrollBody > table > tbody > tr:first-child > td {
                                      border-top: none;
                                    }

                                    div.dataTables_scrollFoot > .dataTables_scrollFootInner {
                                      box-sizing: content-box;
                                    }
                                    div.dataTables_scrollFoot > .dataTables_scrollFootInner > table {
                                      margin-top: 0 !important;
                                      border-top: none;
                                    }

                                    @media screen and (max-width: 767px) {
                                      div.dataTables_wrapper div.dataTables_length,
                                      div.dataTables_wrapper div.dataTables_filter,
                                      div.dataTables_wrapper div.dataTables_info,
                                      div.dataTables_wrapper div.dataTables_paginate {
                                        text-align: center;
                                      }
                                    }
                                    table.dataTable.table-condensed > thead > tr > th {
                                      padding-right: 20px;
                                    }
                                    table.dataTable.table-condensed .sorting:after,
                                    table.dataTable.table-condensed .sorting_asc:after,
                                    table.dataTable.table-condensed .sorting_desc:after {
                                      top: 6px;
                                      right: 6px;
                                    }

                                    table.table-bordered.dataTable th,
                                    table.table-bordered.dataTable td {
                                      border-left-width: 0;
                                    }
                                    table.table-bordered.dataTable th:last-child, table.table-bordered.dataTable th:last-child,
                                    table.table-bordered.dataTable td:last-child,
                                    table.table-bordered.dataTable td:last-child {
                                      border-right-width: 0;
                                    }
                                    table.table-bordered.dataTable tbody th,
                                    table.table-bordered.dataTable tbody td {
                                      border-bottom-width: 0;
                                    }

                                    div.dataTables_scrollHead table.table-bordered {
                                      border-bottom-width: 0;
                                    }

                                    div.table-responsive > div.dataTables_wrapper > div.row {
                                      margin: 0;
                                    }
                                    div.table-responsive > div.dataTables_wrapper > div.row > div[class^="col-"]:first-child {
                                      padding-left: 0;
                                    }
                                    div.table-responsive > div.dataTables_wrapper > div.row > div[class^="col-"]:last-child {
                                      padding-right: 0;
                                    }


                                    @keyframes dtb-spinner {
                                      100% {
                                        transform: rotate(360deg);
                                      }
                                    }
                                    @-o-keyframes dtb-spinner {
                                      100% {
                                        -o-transform: rotate(360deg);
                                        transform: rotate(360deg);
                                      }
                                    }
                                    @-ms-keyframes dtb-spinner {
                                      100% {
                                        -ms-transform: rotate(360deg);
                                        transform: rotate(360deg);
                                      }
                                    }
                                    @-webkit-keyframes dtb-spinner {
                                      100% {
                                        -webkit-transform: rotate(360deg);
                                        transform: rotate(360deg);
                                      }
                                    }
                                    @-moz-keyframes dtb-spinner {
                                      100% {
                                        -moz-transform: rotate(360deg);
                                        transform: rotate(360deg);
                                      }
                                    }
                                    div.dt-button-info {
                                      position: fixed;
                                      top: 50%;
                                      left: 50%;
                                      width: 400px;
                                      margin-top: -100px;
                                      margin-left: -200px;
                                      background-color: white;
                                      border: 2px solid #111;
                                      box-shadow: 3px 3px 8px rgba(0, 0, 0, 0.3);
                                      border-radius: 3px;
                                      text-align: center;
                                      z-index: 21;
                                    }
                                    div.dt-button-info h2 {
                                      padding: 0.5em;
                                      margin: 0;
                                      font-weight: normal;
                                      border-bottom: 1px solid #ddd;
                                      background-color: #f3f3f3;
                                    }
                                    div.dt-button-info > div {
                                      padding: 1em;
                                    }

                                    div.dt-button-collection-title {
                                      text-align: center;
                                      padding: 0.3em 0 0.5em;
                                      font-size: 0.9em;
                                    }

                                    div.dt-button-collection-title:empty {
                                      display: none;
                                    }

                                    ul.dt-button-collection.dropdown-menu {
                                      display: block;
                                      z-index: 2002;
                                      -webkit-column-gap: 8px;
                                      -moz-column-gap: 8px;
                                      -ms-column-gap: 8px;
                                      -o-column-gap: 8px;
                                      column-gap: 8px;
                                    }
                                    ul.dt-button-collection.dropdown-menu.fixed {
                                      position: fixed;
                                      top: 50%;
                                      left: 50%;
                                      margin-left: -75px;
                                      border-radius: 0;
                                    }
                                    ul.dt-button-collection.dropdown-menu.fixed.two-column {
                                      margin-left: -150px;
                                    }
                                    ul.dt-button-collection.dropdown-menu.fixed.three-column {
                                      margin-left: -225px;
                                    }
                                    ul.dt-button-collection.dropdown-menu.fixed.four-column {
                                      margin-left: -300px;
                                    }
                                    ul.dt-button-collection.dropdown-menu > * {
                                      -webkit-column-break-inside: avoid;
                                      break-inside: avoid;
                                    }
                                    ul.dt-button-collection.dropdown-menu.two-column {
                                      width: 300px;
                                      padding-bottom: 1px;
                                      -webkit-column-count: 2;
                                      -moz-column-count: 2;
                                      -ms-column-count: 2;
                                      -o-column-count: 2;
                                      column-count: 2;
                                    }
                                    ul.dt-button-collection.dropdown-menu.three-column {
                                      width: 450px;
                                      padding-bottom: 1px;
                                      -webkit-column-count: 3;
                                      -moz-column-count: 3;
                                      -ms-column-count: 3;
                                      -o-column-count: 3;
                                      column-count: 3;
                                    }
                                    ul.dt-button-collection.dropdown-menu.four-column {
                                      width: 600px;
                                      padding-bottom: 1px;
                                      -webkit-column-count: 4;
                                      -moz-column-count: 4;
                                      -ms-column-count: 4;
                                      -o-column-count: 4;
                                      column-count: 4;
                                    }
                                    ul.dt-button-collection.dropdown-menu .dt-button {
                                      border-radius: 0;
                                    }

                                    div.dt-button-background {
                                      position: fixed;
                                      top: 0;
                                      left: 0;
                                      width: 100%;
                                      height: 100%;
                                      z-index: 2001;
                                    }

                                    @media screen and (max-width: 767px) {
                                      div.dt-buttons {
                                        float: none;
                                        width: 100%;
                                        text-align: center;
                                        margin-bottom: 0.5em;
                                      }
                                      div.dt-buttons a.btn {
                                        float: none;
                                      }
                                    }
                                    div.dt-buttons button.btn.processing,
                                    div.dt-buttons div.btn.processing,
                                    div.dt-buttons a.btn.processing {
                                      color: rgba(0, 0, 0, 0.2);
                                    }
                                    div.dt-buttons button.btn.processing:after,
                                    div.dt-buttons div.btn.processing:after,
                                    div.dt-buttons a.btn.processing:after {
                                      position: absolute;
                                      top: 50%;
                                      left: 50%;
                                      width: 16px;
                                      height: 16px;
                                      margin: -8px 0 0 -8px;
                                      box-sizing: border-box;
                                      display: block;
                                      content: ' ';
                                      border: 2px solid #282828;
                                      border-radius: 50%;
                                      border-left-color: transparent;
                                      border-right-color: transparent;
                                      animation: dtb-spinner 1500ms infinite linear;
                                      -o-animation: dtb-spinner 1500ms infinite linear;
                                      -ms-animation: dtb-spinner 1500ms infinite linear;
                                      -webkit-animation: dtb-spinner 1500ms infinite linear;
                                      -moz-animation: dtb-spinner 1500ms infinite linear;
                                    }
                                </style>'''
                               ]
        self.header_includes = new + list(self.header_includes)
        self.report = report
        # start doc
        self.start_header()
        for line in self.header:
            self.report.write(line)
        self.body = []
        self.update_doc(section_header, center_content, figures, header=True)
    
    def update_doc(self, section_header, center_content,
                   figures, section_footer=False, plotly=False,
                   header=False):
        '''
            Add a new section
        '''
        print('ADDING... ' + section_header)
        self.section_content = []
        self.section_start = []
        #if qc_start:
        #    self.section_start.append('\n# QC Metrics')
        #    self.section_start.append('<details open><summary>Hide QC</summary>')
        if not header:
            self.section_header = '\n# ' + section_header # level 3 headers
            if self.default == 'show':
                self.section_content.append('<details open><summary>Hide</summary>')
                self.section_content.append(center_content)
            elif self.default == 'hide':
                self.section_content.append('<details><summary>Show</summary>')
                self.section_content.append(center_content)
        else:
            self.section_header = ''
            self.section_content = [center_content]
        self.section_start.append(self.section_header)
        if not figures:
            self.figures = []
        else:
            self.figures = figures
        for script, div in self.figures:
            self.section_content.append(div)
            self.section_content.append(script)
        if section_footer:
            self.section_content.append(section_footer)
        if not header:
            self.section_content.append('</details>')
        #if qc_end:
        #    self.section_content.append('</details>')
        self.section_content = self.section_start + self.section_content
        for line in self.add_to_body():
            self.report.write(line)


    def start_header(self):
        '''
            Compose start of header for reports.
            Also prep MD header.
        '''
        with open(self.header_file, 'w') as header_include:
            for line in self.header_includes:
                header_include.write(line + '\n')
        self.header = [line + '\n' for line in self.header]


    def add_to_header(self, header_include_line):
        '''
            Add to header for reports.
        '''
        with open(self.header_file, 'a') as header_include:
            header_include.write(header_include_line + '\n')


    def add_to_body(self):
        '''
            Compose new section.
        '''
        return [section + '\n' for section in self.section_content]


    def write(self):
        '''
            Write out MD file
        '''
        print('Writing...')
        for line in self.header + self.body:
            yield line

    def join_ignore_none(self, lines):
        '''
            join elements of a list (skip None-type objects)
        '''
        return ''.join([line for line in lines if not line == None])






