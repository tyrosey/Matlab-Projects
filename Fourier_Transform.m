%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Tyler Rose
% ECE -- Signals and Systems
%
% Project Description:
% 
% Generate 1D and 2D random signals, calculate the Fourier Transform
% of the generated signals, plot the real and imaginary parts,
% as well as the phases and magnitudes, and then alter the real
% and imaginary parts. Repeat the process for inverse F.T. as well.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% 1-Dimensional signals
close all
 
% Define random vector (1x10)
randomVector_1=rand(1,10).*10;
 
% Define random variable X_0
X_0 = length(randomVector_1);
 
% Equation for Fourier Transform
for k = 1:X_0
    randomVar(k) = 0;
    for n = 1:X_0
        randomVar(k) = randomVar(k) + randomVector_1(n).*exp(((-1j).*2.*pi/X_0).*(n-1).*(k-1));
    end
end

% Creating the plot
figure
subplot(4,1,1)
plot(real(randomVar))
grid on
ylabel('Real part')
subplot(4,1,2)
plot(imag(randomVar))
grid on
ylabel('Imaginary part')
subplot(4,1,3)
plot(abs(randomVar))
grid on
ylabel('Magnitude')
subplot(4,1,4)
grid on
plot(angle_0(randomVar))
ylabel('Phase')
grid on
 
% Double real value, half imaginary value
X_2 = 2.*real(randomVar) + 0.5.*imag(randomVar);
X_3 = length(X_2); 

% Equation for Inverse Fourier Transform
for n=1:X_3
    X_2(n) = 0;
   for k = 1:X_3
      X_2(n) = X_2(n) + ((1/X_3)*(X_2(k)*exp(((1i)*2*pi/X_3)*(k-1)*(n-1)))); 
   end
end
 
figure
stem(real(X_2));
 
% Random Bias of phase
angle_0 = angle_0(X1) + rand(1,length(X1));
magnitude_1 = abs(X1);
X_4 = magnitude_1.*cos(angle_0)+1i.*magnitude_1.*sin(angle_0);
X_5 = length(X_4);
X_6 = zeros(1,length(X_4));
 
% IDTF
for n = 1:X_5
   for k = 1:X_5
      X_6(n) = X_6(n) + ((1/X_5)*(X_4(k)*exp(((1i)*2*pi/X_5)*(k-1)*(n-1)))); 
   end
end
 
figure
stem(real(X_6));
 
 

% 2-Dimensional signals
 
% Define Random Vector (10x10)
randomVector_2=rand(10,10).*10;
X_7 = 10;
 
% Equation for Fourier Transform
for k = 1:10
    X_8(k,k1) = 0;
    for k1 = 1:10
        % rows
        for n = 1:10
            % columns
            for o = 1:10 
                  X_8(k,k1) = X_8(k,k1) + randomVector_2(n,o).*exp(((-1j).*2.*pi *((((n-1)*(k-1))/X_7) + (((o-1)*(k1-1))/X_7))));
            end
        end
    end
end
 
% Creating the plot
figure
subplot(4,2,1)
plot(real(X_8))
grid on
ylabel('Real Value')
subplot(4,2,2)
X4_real_Pcolor = pcolor(real(X_8));
X4_real_Pcolor.FaceColor = 'interp';
title('Pcolor Plot')
grid on
ylabel('Real part')
subplot(4,2,3)
plot(imag(X_8))
grid on
ylabel('Imaginary Value')
subplot(4,2,4)
X4_imag_Pcolor = pcolor(imag(X_8));
X4_imag_Pcolor.FaceColor = 'interp';
title('Pcolor Plot')
grid on
ylabel('Imaginary Value')
subplot(4,2,5)
plot(abs(X_8))
grid on
ylabel('Magnitude')
subplot(4,2,6)
X4_abs_Pcolor = pcolor(abs(X_8));
X4_abs_Pcolor.FaceColor = 'interp';
title('Pcolor Plot')
grid on
ylabel('Magnitude')
subplot(4,2,7)
grid on
plot(angle_0(X_8))
ylabel('Phase')
subplot(4,2,8)
grid on
X4_angle_Pcolor = pcolor(angle_0(X_8));
X4_angle_Pcolor.FaceColor = 'interp';
title('Pcolor Plot')
grid on
ylabel('Phase')
grid on
 
% Double real value, half imaginary value
X_9 = 2 .* real(X_8) + 0.5 .* imag(X_8);
X_10 = 10;
 
% Equation for Inverse Fourier Transform
for n = 1:10
    for o = 1:10
        X_11(n,o) = 0;
        for k = 1:10
            for k1 = 1:10
                  X_11(n,o) = X_11(n,o) + (1 ./ X_10) *(X_9(k,k1).*exp(((1j).*2.*pi *((((n-1)*(k-1))/X_10)+ (((o-1)*(k1-1))/X_10)))));
            end
        end
    end
end
 
figure(5);
hold on;
subplot(1,2,1)
stem(real(X_11));
title('INVERSE DFT');
subplot(1,2,2)
x5_pcolor = pcolor(real(X_11));
x5_pcolor.FaceColor = 'interp';
title('Pcolor Plot')
hold off;
 
% Random Bias of phase
angle_2 = angle_0(X_8)+rand(1,length(X_8));
magnitude_2 = abs(X_8);
X_12 = magnitude_2.*cos(angle_2) + 1j.*magnitude_2.*sin(angle_2);
X_13 = 10;
 
 
% IDTF
for n = 1:10
    for o = 1:10
        X_14(n,o) = 0;
        for k = 1:10
            for k1 = 1:10
                  X_14(n,o) = X_14(n,o) + (1 ./ X_13) *(X_12(k,k1).*exp(((1j).*2.*pi *((((n-1)*(k-1))/X_13)+ (((o-1)*(k1-1))/X_13)))));
            end
        end
    end
end
 
% Plot
figure(6);
hold on;
subplot(1,2,1);
stem(real(X_14));
title('IDFT');
subplot(1,2,2)
x6_pcolor = pcolor(real(X_14));
x6_pcolor.FaceColor = 'interp';
title('Pcolor Plot')
 
hold off;
