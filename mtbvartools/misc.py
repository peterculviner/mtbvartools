import subprocess, json
import pandas as pd
import numpy as np

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