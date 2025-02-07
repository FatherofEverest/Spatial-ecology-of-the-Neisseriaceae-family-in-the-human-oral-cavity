#!/bin/bash
METAPHLAN_4_DIR=/workspace/jmarkwelchlab/P_0622_Haemophilus_Aggregatibacter/METAPHLAN_4
sample2markers.py -i $METAPHLAN_4_DIR/sams/*.sam.bz2 \
  -o $METAPHLAN_4_DIR/consensus_markers \
  --tmp $METAPHLAN_4_DIR/sample2markersTMP \
  --debug \
  -n 20

