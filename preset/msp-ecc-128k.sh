#!/bin/bash
CC="msp430-gcc -mmcu=msp430f1611" CXX="c++" cmake -DARITH=msp-asm -DCMAKE_SYSTEM_NAME=Generic -DALIGN=2 -DARCH=MSP -DBENCH=1 "-DBN_METHD=BASIC;MULTP;MONTY;BASIC;BASIC;BASIC" -DCHECK=OFF -DCOLOR=OFF "-DCFLAGS:STRING=-O2 -g -mmcu=msp430f1611 -ffunction-sections -fdata-sections -fno-inline -mdisable-watchdog" -DDOCUM=OFF -DEB_DEPTH=3 -DEB_PLAIN=OFF -DEP_DEPTH=3 -DEP_PLAIN=OFF "-DFB_METHD=LODAH;RLC_TABLE;QUICK;QUICK;BASIC;BASIC;EXGCD;BASIC;BASIC" -DFB_PRECO=OFF "-DFP_METHD=BASIC;COMBA;COMBA;QUICK;LOWER;BASIC" -DFP_PMERS=ON "-DLDFLAGS=-Wl,--gc-sections" -DSEED= -DSHLIB=OFF -DSTRIP=ON -DTESTS=1 -DTIMER=CYCLE -DVERBS=OFF -DWSIZE=16 -DFP_PRIME=256 -DFB_POLYN=283 -DBN_PRECI=284 -DMD_METHD=SH256 "-DWITH=FP;FB;EP;EB;EC;DV;CP;MD;BN" -DEC_METHD=PRIME -DRAND=HASHD $1
