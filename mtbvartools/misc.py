import subprocess, timeit, json, os, pickle, re
import pandas as pd
import numpy as np
from datetime import datetime, timedelta

def contShell(cmd, is_return=False, exception_at=None):
    def yieldcmd(cmd):
        popen = subprocess.Popen(
            cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        for stdout_line in iter(popen.stdout.readline, ""):
            yield stdout_line
        popen.stdout.close()
        return_code = popen.wait()
        if return_code:
            raise subprocess.CalledProcessError(return_code, cmd)
    if not is_return:
        for line in yieldcmd(cmd):
            if exception_at is not None and exception_at in line[:-1]:
                raise ValueError(f'Found "{exception_at}" in output\n {line[:-1]}.')
            print(line[:-1])
    if is_return:
        output = [line for line in yieldcmd(cmd)]
        if len(output) > 0:
            return '\n'.join(output)
        else:
            return ''


def getSlices(step, length):
    slice_intervals = []
    left = 0
    while left < length:
        slice_intervals.append((left, left + step))
        left += step
    return slice_intervals


def chainLookup(keys, dictionary):
        try:
            value = dictionary
            for key in keys:
                value = value[key]
            return value
        except KeyError:
            return np.nan


def getJSONValues(filename, keys_list):
    try:
        loaded_json = json.load(open(filename, 'r'))
    except FileNotFoundError:
        return pd.Series(
            {label: np.nan for label, keys in keys_list})
    return pd.Series(
        {label: chainLookup(keys, loaded_json) for label, keys in keys_list})

def getNearbyGenes(pos, gene_table):
    output_dict = {}
    # ZERO INDEXED POS
    # prepare values for inside genes
    inside = gene_table.index[np.all([gene_table.start <= pos, gene_table.end >= pos], axis=0)]
    dinside = []
    for geneid in inside:
        if gene_table.loc[geneid, 'strand'] == 0:
            dinside.append(str(gene_table.loc[geneid, 'start'] - pos))
        elif gene_table.loc[geneid, 'strand'] == 1:
            dinside.append(str(gene_table.loc[geneid, 'end'] - pos))
    output_dict['inside'] = ','.join(inside)
    output_dict['dinside'] = ','.join(dinside)
    # left gene
    try:
        output_dict['left'] = gene_table.index[gene_table.end < pos][-1]
        output_dict['dleft'] = str(gene_table.loc[output_dict['left'], 'end'] - pos)
    except IndexError:
        output_dict['left'] = ''
        output_dict['dleft'] = ''
    # right gene
    try:
        output_dict['right'] = gene_table.index[gene_table.start > pos][0]
        output_dict['dright'] = str(gene_table.loc[output_dict['right'], 'start'] - pos)
    except IndexError:
        output_dict['right'] = ''
        output_dict['dright'] = ''
    return output_dict

class StopWatch():
    def start(self, name, flat_file=None):
        if name not in self.ignore:
            now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            start = timeit.default_timer()
            self.timestamps[name] = {'start': now, 'tstart': start}
            out_str = f'{now} - START {name}'
            print(out_str)
            if flat_file is not None:
                with open(flat_file, 'a') as f:
                    f.write(out_str + '\n')
            if self.dict_path is not None:
                with open(self.dict_path, 'wb') as f:
                    pickle.dump(self.timestamps, f)
         
    def end(self, name, flat_file=None, notes=None):
        if name not in self.ignore:
            if name not in self.timestamps.keys():
                raise KeyError(f'{name} not found in tracked timestamps.')
            now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
            end = timeit.default_timer()
            record = self.timestamps[name]
            record['end'] = now
            record['tend'] = end
            record['selapsed'] = round(timeit.default_timer() - record['tstart'])
            record['elapsed'] = str(timedelta(seconds=record['selapsed']))
            out_str = f'{now} - END {name} elapsed {record["selapsed"]} ({record["elapsed"]})'
            if notes is not None:
                record['notes'] = notes
                out_str += f' NOTE: {notes}'
            print('\n' + out_str + '\n')
            if flat_file is not None:
                with open(flat_file, 'a') as f:
                    f.write(out_str + '\n')
            if self.dict_path is not None:
                with open(self.dict_path, 'wb') as f:
                    pickle.dump(self.timestamps, f)

    def report(self, file):
        with open(file, 'w') as f:
            f.write(f'# starts={self.timestamps["info"]["starts"]}\n')
            f.write(f'# cpu={",".join(self.timestamps["info"]["cpu"])}\n')
            pd.DataFrame(
                data={k: v for k, v in self.timestamps.items() if k != 'info'}).T.to_csv(f)
        
    def __init__(self, dict_path=None):
        self.dict_path = dict_path
        if dict_path is not None and os.path.exists(dict_path):
            with open(dict_path, 'rb') as f:
                self.timestamps = pickle.load(f)
            self.timestamps['info']['starts'] += 1
            # already started once, do not overwrite completed timings
            self.ignore = []
            for key, value in self.timestamps.items():
                if 'end' in value.keys():  # already complete, ignore
                    self.ignore.append(key)
        else:
            self.timestamps = {'info':{'starts': 1}}
            self.timestamps['info']['cpu'] = []
            self.ignore = []
        # store hardware information
        self.timestamps['info']['cpu'].append(
            re.sub('Model name:\s*', '', subprocess.check_output("lscpu | grep 'Model name:'", shell=True).strip().decode()))