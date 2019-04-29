// Name of Experiment: Experiment #20190409023153 v10

OPENQASM 2.0;
include "qelib1.inc";


qreg q[5];
creg c[5];

h q[0];
h q[1];
cx q[1],q[0];
x q[0];
h q[2];
cx q[2],q[0];
h q[0];
cx q[1],q[0];
x q[0];
cx q[2],q[0];
measure q[2] -> c[2];
measure q[1] -> c[1];
