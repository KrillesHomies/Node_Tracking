function [fg, fb, fd, tg, tb ,td] = calculate_forces(state, shape, ocean, Constants,step_rv, T)

% Gravity 
fg = Constants.g'*shape.Mass;
tg = cross(q2dcm(state.Q)'*shape.Center_of_mass,fg);

% Bouyancy
if(step_rv == 0)
    fb = Constants.p*(-Constants.g')*Constants.v_sinking;
else
    fb = Constants.p*(-Constants.g')*Constants.v_rising;
end
% Incase node is above water
if(state.P(end) < 0)
    fb = fb/4;
end
if(step_rv == 0)
    tb = cross(q2dcm(state.Q)'*shape.Center_of_Bouyancy_sink,fb);
else
    tb = cross(q2dcm(state.Q)'*shape.Center_of_Bouyancy_rise,fb);
end


% FVelocity from current in the water
%Rotation 
rot = ocean.current_rotation_intial + [0; 0; sin(2*pi*T/ocean.rotation_cycle_time)] + normrnd([0 0 0]',ocean.current_rotation_noise);
% Amplitute
V_c = ocean.current_amplitude_intial + ocean.current_amplitude_timevarying*sin(2*pi*T/ocean.amplitude_cycle_time) + normrnd([0 0 0]',ocean.current_amplitude_noise);
V_c =  eul2rotm(rot','XYZ')*V_c;

% Drag forces 
if(step_rv == 0)
    fd = 0.5*Constants.p*Constants.C_sinking*(state.V - V_c).^2.*sign(V_c - state.V);
else
    fd = 0.5*Constants.p*Constants.C_rising*(state.V - V_c).^2.*sign(V_c - state.V);
end
td = cross(q2dcm(state.Q)'*shape.Center_of_Pressure,fd);