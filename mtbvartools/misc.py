import subprocess, json, timeit
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

class StopWatch():
    def start(self, name, file=None):
        now = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        start = timeit.default_timer()
        self.timestamps[name] = {'start': now, 'tstart': start}
        out_str = f'{now} - START {name}'
        print(out_str)
        if file is not None:
            with open(file, 'a') as f:
                f.write(out_str + '\n')
         
    def end(self, name, file=None, notes=None):
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
        if file is not None:
            with open(file, 'a') as f:
                f.write(out_str + '\n')

    def report(self, file):
        pd.DataFrame(data=self.timestamps).T.to_csv(file)
        
    def __init__(self):
        self.timestamps = {}