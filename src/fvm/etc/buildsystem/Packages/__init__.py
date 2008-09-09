def _listModules():
    import os
    mods = []
    for filename in os.listdir(__path__[0]):
        base, ext = os.path.splitext(filename)
        if ext == '.py' and base != '__init__':
            mods.append(base)
    return mods

exec('import %s' % ','.join(_listModules()))
