center_forces = [0.427910390000000,0.605873072000000,0.375153606000000,0.201349709000000,0.0964100200000000;2.33197180300000,0.671104802000000,0.463204410000000,0.204180681000000,0.0496439310000000;0.319765419000000,0.161184061000000,0.0792299040000000,0.00331016900000000,0.00443828400000000;0.122823840000000,0.0348303060000000,0.00526388900000000,0.00392292600000000,0.00274043800000000;0.0649549410000000,0.00303949600000000,0.00274043800000000,0.00235375600000000,0.00235375600000000];
frequencies = [61,434,976.600000000000,1953.10000000000,7812.50000000000];
duty_cycles = [50,10,5,1,0.500000000000000];

center_forces = fliplr(center_forces);
duty_cycles = fliplr(duty_cycles);

surf(duty_cycles,frequencies,center_forces)
xlabel('Duty Cycles (%)')
ylabel('Frequencies (Hz)')
zlabel('Force (g)')

