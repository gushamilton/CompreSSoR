import json, sys, os
root, out = sys.argv[1:]
idx = json.load(open(os.path.join(root, 'native.index.json')))
# These are the actual physical stream types in the 0.4.4 native store:
# position=u32, substitution/eaf/se=u8, z=u16.  Do not infer this from the
# logical semantic bit widths in manifest.json.
types = {'position':1, 'substitution':10, 'z':7, 'eaf':10, 'se':10}
with open(out, 'w') as f:
    f.write('stream\tfile\tdtype\toffset\tlength\tvalues\trow_start\trow_stop\tfirst_position\tlast_position\tblock_id\n')
    for stream, spec in idx['streams'].items():
        blocks = spec['blocks']
        for j,b in enumerate(blocks):
            if stream in ('position','substitution'):
                # Map 8,192-row key blocks to their enclosing 65,536-row value block.
                block_id = b['row_start'] // idx['value_block_rows']
            else:
                block_id = j
            # Index has matching genomic bounds for both key and value block families.
            if stream in ('position','substitution'):
                kb = idx['key_blocks'][j]
                first,last = kb['first_position'],kb['last_position']
            else:
                vb = idx['value_blocks'][j]
                first,last = vb['first_position'],vb['last_position']
            path=os.path.join(root,spec['file'])
            f.write(f"{stream}\t{path}\t{types[stream]}\t{b['offset']}\t{b['length']}\t{b['values']}\t{b['row_start']}\t{b['row_stop']}\t{first}\t{last}\t{block_id}\n")
