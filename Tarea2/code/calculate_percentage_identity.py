import sys
def amounts_of_shared_structures(known, predict, show_percentage = True):
    assert len(known) == len(predict)
    
    shared_amounts = {'E': 0, 'C': 0, 'H': 0}
    for structure1, structure2 in zip(known, predict):

        if structure1 == structure2:
            try:
                shared_amounts[structure1] += 1
            except KeyError as error:
                if structure1 != '-':
                    raise ValueError(error)
    
    if show_percentage:
        # Side Effect.
        total = sum(shared_amounts.values()) * 100 / len(known)
        print("Global identity percentage: {}%".format(total))
        for key, value in shared_amounts.items():
            try:
                print('Itentity percentage of {}: {}%'.format(key, value * 100 / predict.count(key)))
            except ZeroDivisionError:
                print("Identity percentage of {}: 0.0%".format(key))
    
    return(shared_amounts)

  
# Run function with our sequences.
sequences_file = sys.argv[1]
with open(sequences_file, "r") as f:
  lines = f.readlines()
  known, predict = lines[3].strip('\n'), lines[5].strip('\n')

amounts_of_shared_structures(known, predict, show_percentage = True)