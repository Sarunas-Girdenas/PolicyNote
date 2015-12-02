% in this file we compute IRFs figures

% compute the first figure with IRFs

load figure1_data

figure(1)
subplot(2,2,1)
plot(figure1_data.u_br)
hold
plot(figure1_data.u_u,'k')
plot(figure1_data.u_o,'r')
title('Unemployment')
xlabel('Time Periods')
subplot(2,2,2)
plot(figure1_data.i_br)
hold
plot(figure1_data.i_u,'k')
plot(figure1_data.i_o,'r')
title('Interest Rate')
xlabel('Time Periods')
subplot(2,2,3)
plot(figure1_data.b_br)
hold
plot(figure1_data.b_u,'k')
plot(figure1_data.b_o,'r')
title('Borowing')
xlabel('Time Periods')
subplot(2,2,4)
plot(figure1_data.infl_br)
hold
plot(figure1_data.infl_u,'k')
plot(figure1_data.infl_o,'r')
title('Inflation')
xlabel('Time Periods')

clear

load figure2_data

figure(2)
subplot(2,2,1)
plot(figure2_data.w_br)
hold
plot(figure2_data.w_u,'k')
plot(figure2_data.w_o,'r')
title('Wage')
xlabel('Time Periods')
subplot(2,2,2)
plot(figure2_data.Y_br)
hold
plot(figure2_data.Y_u,'k')
plot(figure2_data.Y_o,'r')
title('Output')
xlabel('Time Periods')
subplot(2,2,3)
plot(figure2_data.S_br)
hold
plot(figure2_data.S_u,'k')
plot(figure2_data.S_o,'r')
title('Surplus')
xlabel('Time Periods')
subplot(2,2,4)
plot(figure2_data.q_br)
hold
plot(figure2_data.q_u,'k')
plot(figure2_data.q_o,'r')
title('VF Probability')
xlabel('Time Periods')

