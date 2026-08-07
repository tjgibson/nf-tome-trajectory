#!/usr/bin/python3

import json
import sys


time_point = sys.argv[1].split(",")
tree_edge_fn = sys.argv[2]
edge_prob_fn = sys.argv[3]
celltype_groups_fn = sys.argv[4]


time_point_n = len(time_point)
time_point_id = {}
for i in range(0,time_point_n):
    time_point_id[time_point[i]] = i


### read edge info 
edge = {}; node_all = []; node_each = {}; main_edge = set()

file = open(tree_edge_fn)
for line in file:
    l = line.rstrip().split('\t')

    if l[0] not in node_all:
        node_all.append(l[0])
    if l[1] not in node_all:
        node_all.append(l[1])

    edge[l[0]] = edge.get(l[0], [])
    edge[l[0]].append(l[1])

    num = time_point_id[l[0].split(':')[0]]
    node_each[num] = node_each.get(num, [])
    if l[0] not in node_each[num]:
        node_each[num].append(l[0])

    main_edge.add((l[0],l[1]))

### read edge weight file
main_weight = {}; extra_weight = {}
file = open(edge_prob_fn)
for line in file:
    l = line.rstrip().split('\t')
    if float(l[2])>0.8:
        xx = "specific"
    elif float(l[2])>0.2:
        xx = "additional"
    if (l[0],l[1]) in main_edge:
        if float(l[2]) > 0.2:
            main_weight[l[1]] = float(l[2])
        else:
            main_weight[l[1]] = 0
    else:
        if float(l[2]) > 0.2:
            extra_weight[(l[0],l[1])] = float(l[2])


### read group information, used for the color of the node
file = open(celltype_groups_fn)
i = 1
coor = {}
node_group = {}
color_map = {}
for line in file:
    l = line.rstrip().split("\t")
    node_group[l[0]] = int(l[1])
    coor[l[0]] = i
    color_map[int(l[1])] = l[2]
    i += 1
file.close()

### create the info for all the node
dat = {}

for i in node_all:

    dat[i] = {'name':i}
    if i in edge:
        dat[i]['children'] = []

    if i in main_weight:
        dat[i]['edge_weight'] = main_weight[i]
    else:
        print(i)

    ### 10. coors of x (used to create the plot, the cell type order)
    if i.split(':')[1] in coor:
        dat[i]['fx'] = str(coor[i.split(':')[1]])
    else:
        print(i)


    if i.split(':')[1] in node_group:
        dat[i]['node_group'] = color_map[node_group[i.split(':')[1]]]
    else:
        print(i)


### add extra node first
for i in extra_weight:
    tmp = {'name': i[1], 'edge_weight': extra_weight[i], 'fx': str(coor[i[1].split(':')[1]]), 'node_group': color_map[node_group[i[1].split(':')[1]]]}
    dat[i[0]]['children'] = dat[i[0]].get('children',[])
    dat[i[0]]['children'].append(tmp)
    print(dat[i[0]]['children'])



### connect each time point
for i in range(time_point_n-2, 0, -1):
    for j in node_each[i]:
        for k in edge[j]:
            dat[j]['children'].append(dat[k])

for k in edge['E3:Morula']:
    dat['E3:Morula']['children'].append(dat[k])
dat_json = dat['E3:Morula']


with open("./tree.json", 'w') as json_file:
    json.dump(dat_json, json_file)