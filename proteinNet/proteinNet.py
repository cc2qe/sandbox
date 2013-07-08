#!/usr/bin/env python

import sys, argparse

def proteinNet(proteinA, proteinB, proteinFile):
    # load the proteinFile into an array
    protein_pairs = []
    protein_set = set()
    protein_hash = dict()
    for line in proteinFile:
        v = line.rstrip().split('\t')
        if not v[0] in protein_hash:
            protein_hash[v[0]] = list()
        if not v[1] in protein_hash:
            protein_hash[v[1]] = list()

        # if protein is not paired with itself, add to the list of binary interactions
        if v[0] != v[1]:
            protein_hash[v[0]].append(v[1])
            protein_hash[v[1]].append(v[0])

#            protein_pairs.append([v[0], v[1]])
            protein_set.add(v[0])
            protein_set.add(v[1])        

    if not (proteinA in protein_set and proteinB in protein_set):
        sys.stderr.write('One or both proteins absent from database\n')
        sys.exit(0)
    
    # chain_list is a list of chains that we'll try to get from proteinA to proteinB
    chain_list = [[proteinA]]

    # the genes that have already been visited
    visited = set()
    while 1:
        # protein interaction partner from this round of iterations
        new_partners = []
        to_delete = list()
        for chain in chain_list:
            markForDeletion = True
            
            last_in_chain = chain[-1]

            if last_in_chain in protein_hash:
                # the next link in the chain
                next_protein_list = protein_hash[last_in_chain]

                if proteinB in next_protein_list:
                    chain.append(proteinB)
                    print '\t'.join(chain)
                    sys.exit(0)
                else:
                    for next_prot in next_protein_list:
                        if not next_prot == last_in_chain and not next_prot in chain and not next_prot in visited:
                            new_partners.append(next_prot)
                            visited.add(next_prot)
                            markForDeletion = False
                    if markForDeletion:
                        to_delete.append(chain)

        for del_chain in to_delete:
            chain_list.remove(del_chain)
            # print 'removed', del_chain

        new_chain_list = []
        for chain in chain_list:
            for partner in new_partners:
                temp_chain = list(chain)
                temp_chain.append(partner)
                if partner == proteinB:
                    print temp_chain, len(temp_chain)
                    sys.exit(0)
                new_chain_list.append(temp_chain)

        # save the new chain list as chain_list
        chain_list = new_chain_list
        # print len(chain_list)


# --------------------------------------
# argument parsing

def main():
    parser = argparse.ArgumentParser(description='get distance between set of n genes based on protein protein interaction')
    parser.add_argument('-a', '--proteinA', type=str, required=False)
    parser.add_argument('-b', '--proteinB', type=str, required=False)
    parser.add_argument('-l', '--proteinList', type=str, required=False, help='comma separated list of proteins')
    parser.add_argument('-p', '--proteinFile', type=argparse.FileType('r'), required=True, help='Binary protein interaction file (\'-\' for stdin)')

    args = parser.parse_args()

    if args.proteinList != None:
        pList = args.proteinList.rstrip().split(',')
        (args.proteinA, args.proteinB) = pList

    proteinNet(args.proteinA, args.proteinB, args.proteinFile)

if __name__ == '__main__':
    sys.exit(main())

