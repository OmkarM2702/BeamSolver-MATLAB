clc;
close all;

disp("Beam Analysis Program");
disp("1. Point Load");
disp("2. Uniformly Distributed Load (UDL)");
disp("3. Uniformly Varying Load (UVL)");
disp("4. Moment");
choice = input("Enter your choice (1-4): ");

switch choice
    case 1
        % Point Load
        len = input("Enter length of the beam (m): ");
        disp("Enter '1' for cantilever, '2' for simply supported");
        type = input("Beam type: ");
        mag = input("Enter magnitude of point load (N): ");
        loc = input("Enter location of load from left (m): ");

        if type == 1
            % Cantilever beam
            Ra = mag;
            M = mag * loc;

            % SFD
            figure;
            x = [0 loc len];
            V = [-Ra -Ra 0];
            plot(x, V, 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Cantilever with Point Load');

            % BMD
            figure;
            x = [0 loc len];
            M = [0 -mag*loc -mag*loc];
            plot(x, M, 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Cantilever with Point Load');

        elseif type == 2
            % Simply supported beam
            Ra = mag * (len - loc) / len;
            Rb = mag - Ra;

            % SFD
            figure;
            x = [0 loc len];
            V = [Ra Ra -Rb];
            plot(x, V, 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Simply Supported with Point Load');

            % BMD
            figure;
            x = 0:0.01:len;
            M = Ra*x;
            M(x >= loc) = Ra*loc - mag*(x(x >= loc) - loc);
            plot(x, M, 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Simply Supported with Point Load');
        end

    case 2
        % Uniformly Distributed Load (UDL)
        len = input("Enter length of beam (m): ");
        disp("Enter '1' for cantilever, '2' for simply supported");
        type = input("Beam type: ");
        mag = input("Enter intensity of UDL (N/m): ");
        loc1 = input("Enter start of UDL (m): ");
        loc2 = input("Enter end of UDL (m): ");

        w = mag;
        L = len;

        if type == 1
            % Cantilever
            total_load = w * (loc2 - loc1);
            eq_point = loc1 + (loc2 - loc1)/2;
            M_eq = total_load * eq_point;

            % SFD
            figure;
            x1 = 0:0.01:loc1;
            sf1 = zeros(1, length(x1));
            x2 = loc1:0.01:loc2;
            sf2 = -w * (x2 - loc1);
            x3 = loc2:0.01:len;
            sf3 = -total_load * ones(1, length(x3));
            plot([x1 x2 x3], [sf1 sf2 sf3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Cantilever with UDL');

            % BMD
            figure;
            b1 = 0:0.01:loc1;
            bm1 = zeros(1, length(b1));
            b2 = loc1:0.01:loc2;
            bm2 = -w*(b2 - loc1).^2 / 2;
            b3 = loc2:0.01:len;
            bm3 = -M_eq * ones(1, length(b3));
            plot([b1 b2 b3], [bm1 bm2 bm3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Cantilever with UDL');

        elseif type == 2
            % Simply Supported
            total_load = w * (loc2 - loc1);
            x_bar = loc1 + (loc2 - loc1)/2;
            Rb = -(total_load * x_bar) / len;
            Ra = total_load + Rb;

            % SFD
            figure;
            x1 = 0:0.01:loc1;
            sf1 = Ra * ones(1, length(x1));
            x2 = loc1:0.01:loc2;
            sf2 = Ra - w * (x2 - loc1);
            x3 = loc2:0.01:len;
            sf3 = (Ra - total_load) * ones(1, length(x3));
            plot([x1 x2 x3], [sf1 sf2 sf3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Simply Supported with UDL');

            % BMD
            figure;
            b1 = 0:0.01:loc1;
            bm1 = Ra * b1;
            b2 = loc1:0.01:loc2;
            bm2 = Ra * b2 - w*(b2 - loc1).^2 / 2;
            b3 = loc2:0.01:len;
            bm3 = Ra * b3 - total_load * (b3 - x_bar);
            plot([b1 b2 b3], [bm1 bm2 bm3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Simply Supported with UDL');
        end

    case 3
        % Uniformly Varying Load (UVL)
        len = input("Enter length of beam (m): ");
        disp("Enter '1' for cantilever, '2' for simply supported");
        type = input("Beam type: ");
        w1 = input("Enter load at start (N/m): ");
        w2 = input("Enter load at end (N/m): ");
        loc1 = input("Enter start location of UVL (m): ");
        loc2 = input("Enter end location of UVL (m): ");

        if type == 1
            % Cantilever
            x = loc1:0.01:loc2;
            wx = w1 + ((w2 - w1)/(loc2 - loc1)) * (x - loc1);
            area = trapz(x, wx);
            centroid = trapz(x, x .* wx) / area;

            % SFD
            figure;
            x1 = 0:0.01:loc1;
            sf1 = zeros(1, length(x1));
            x2 = loc1:0.01:loc2;
            wx2 = w1 + ((w2 - w1)/(loc2 - loc1)) * (x2 - loc1);
            sf2 = -cumtrapz(x2, wx2);
            x3 = loc2:0.01:len;
            sf3 = sf2(end) * ones(1, length(x3));
            plot([x1 x2 x3], [sf1 sf2 sf3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Cantilever with UVL');

            % BMD
            figure;
            b1 = 0:0.01:loc1;
            bm1 = zeros(1, length(b1));
            b2 = loc1:0.01:loc2;
            bm2 = -cumtrapz(b2, cumtrapz(b2, wx2));
            b3 = loc2:0.01:len;
            bm3 = bm2(end) * ones(1, length(b3));
            plot([b1 b2 b3], [bm1 bm2 bm3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Cantilever with UVL');

        elseif type == 2
            % Simply Supported
            x = loc1:0.01:loc2;
            wx = w1 + ((w2 - w1)/(loc2 - loc1)) * (x - loc1);
            area = trapz(x, wx);
            x_bar = trapz(x, x .* wx) / area;
            Rb = -area * x_bar / len;
            Ra = area + Rb;

            % SFD
            figure;
            x1 = 0:0.01:loc1;
            sf1 = Ra * ones(1, length(x1));
            x2 = loc1:0.01:loc2;
            wx2 = w1 + ((w2 - w1)/(loc2 - loc1)) * (x2 - loc1);
            sf2 = Ra - cumtrapz(x2, wx2);
            x3 = loc2:0.01:len;
            sf3 = sf2(end) * ones(1, length(x3));
            plot([x1 x2 x3], [sf1 sf2 sf3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Simply Supported with UVL');

            % BMD
            figure;
            b1 = 0:0.01:loc1;
            bm1 = Ra * b1;
            b2 = loc1:0.01:loc2;
            bm2 = Ra * b2 - cumtrapz(b2, cumtrapz(b2, wx2));
            b3 = loc2:0.01:len;
            bm3 = bm2(end) + (Ra - area) * (b3 - loc2);
            plot([b1 b2 b3], [bm1 bm2 bm3], 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Simply Supported with UVL');
        end

    case 4
        % Moment
        len = input("Enter length of beam (m): ");
        disp("Enter '1' for cantilever, '2' for simply supported");
        type = input("Beam type: ");
        M = input("Enter applied moment (N-m): ");
        loc = input("Enter location of moment from left (m): ");

        if type == 1
            % Cantilever
            figure;
            x = [0 loc len];
            SF = [0 0 0];
            plot(x, SF, 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Cantilever with Moment');

            figure;
            BM = [0 -M -M];
            plot(x, BM, 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Cantilever with Moment');

        elseif type == 2
            % Simply Supported
            Ra = M / len;
            Rb = -Ra;

            figure;
            x = [0 loc len];
            SF = [Ra Ra Rb];
            plot(x, SF, 'm'); grid on;
            xlabel('x (m)'); ylabel('Shear Force (N)');
            title('SFD - Simply Supported with Moment');

            figure;
            x = 0:0.01:len;
            BM = Ra * x;
            BM(x >= loc) = Ra * x(x >= loc) - M;
            plot(x, BM, 'm'); grid on;
            xlabel('x (m)'); ylabel('Bending Moment (N-m)');
            title('BMD - Simply Supported with Moment');
        end

    otherwise
        disp("Invalid choice.");
end
