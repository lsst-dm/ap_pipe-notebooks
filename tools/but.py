def list_datasets(butler,collections,dataID,datasetType,where='',max_print=20,
                  print_run=True,print_type=False,print_ID=True):
    
    results = butler.registry.queryDatasets(collections=collections,
                                            dataId=dataID,
                                            datasetType=datasetType,
                                            where=where,
                                           )

    results = list(results)
    L = len(results)
    print(f'Found {L}')
    for i in range(min(L,max_print)):
        r = results[i]
        
        s = ''
        if print_run:
            s += r.run + ' '
        if print_type:
            s += f'{r.datasetType.name}' + ' '
        if print_ID:
            s += f'{r.dataId.full}'
        print(s)

    return results
