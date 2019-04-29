// Name of Experiment: Experiment #20190319143223 v1

OPENQASM 2.0;
include "qelib1.inc";


qreg q[5];
creg c[5];

h q[1];
h q[2];
h q[0];
cx q[2],q[0];
cx q[1],q[0];
cx q[2],q[0];
h q[0];
h q[2];
h q[0];
h q[1];
x q[0];
x q[1];
h q[0];
cx q[1],q[0];
h q[0];
cx q[1],q[0];
h q[0];
x q[0];
x q[1];
h q[0];
h q[1];
measure q[2] -> c[2];
measure q[1] -> c[1];
measure q[0] -> c[0];
