"""JSON-safe solver diagnostics, separate from geometry and numeric inputs."""
import math
import numbers


def json_safe_solver_diagnostics(receipt):
    """Retain exact status and mark non-finite diagnostic values explicitly.

    This must only be used on solver *reports*. It does not accept or repair
    geometry, constraints, a candidate, or a solution vector. Strict JSON
    serialization remains appropriate for those numerical inputs and outputs.
    """
    if not isinstance(receipt, dict):
        raise TypeError('A solver diagnostic dictionary is required')
    if not isinstance(receipt.get('status'), str):
        raise ValueError('The exact string solver status must be retained')
    reserved={'nonfiniteDiagnosticValues', 'nonfiniteDiagnosticFields'}
    if reserved.intersection(receipt):
        raise ValueError('Solver receipt already contains sanitization metadata')
    nonfinite=[]

    def clean(value, path):
        if value is None or isinstance(value, (str, bool)):
            return value
        if isinstance(value, numbers.Integral):
            return int(value)
        if isinstance(value, numbers.Real):
            number=float(value)
            if math.isfinite(number):
                return number
            classification='nan' if math.isnan(number) else ('+infinity' if number>0 else '-infinity')
            nonfinite.append({'path':path,'classification':classification})
            return None
        if isinstance(value, dict):
            if not all(isinstance(key,str) for key in value):
                raise TypeError('Solver diagnostic keys must be strings')
            return {key:clean(item,path+'/'+key.replace('~','~0').replace('/','~1')) for key,item in value.items()}
        if isinstance(value, (list, tuple)):
            return [clean(item,path+'/'+str(index)) for index,item in enumerate(value)]
        raise TypeError('Unsupported solver diagnostic type: '+type(value).__name__)

    result=clean(receipt,'')
    assert result['status']==receipt['status']
    result['nonfiniteDiagnosticValues']=bool(nonfinite)
    result['nonfiniteDiagnosticFields']=nonfinite
    return result
