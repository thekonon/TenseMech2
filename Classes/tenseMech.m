classdef tenseMech<TensegritySettings
    %Veřejné metody
    properties(Access = public)
        Y_sim
        t_sim
    end
    %Privátní metody
    properties(Access = private)
        M
        Phi
        PhiD
        W
        Tc
        cable_forces
        residuum
        megaMatrixFK
        megaRightSideFK

        nodes_velocity


        %Stavové proměnné
        s
        sd

    end
    %Konstanty
    properties(Access = public,Constant)
        time_stop = 5
        alf = 1
        bet = 1
        g = -9.81*0
    end

    methods(Access = public)
        function obj = tenseMech()
            obj@TensegritySettings()
            obj.initialConditionsToStates()
            obj.stateToNodes()
            obj.solveFK()
        end
        function solveFK(obj, timeStop)
            if nargin == 1
                timeStop = obj.time_stop;
            end
            options = odeset('RelTol',1e-4,'AbsTol',1e-5);
            [t, Y] = ode15s(@obj.stepFK, [0, timeStop], [obj.s; obj.sd],options);
            obj.Y_sim = Y;
            obj.t_sim = t;
            figure(1), clf,
            for i = 1:1:numel(t)-1
                obj.s = Y(i, 1:42)';
                obj.stateToNodes()
                obj.plotNodes
                xlim([-0.200 0.200])
                ylim([-0.200 0.200])
                zlim([-0.600 0.600])
                view([-180.900 21.200]) %normální pohled
%                 view([-180.049 90.000]) %pohled dolu
%                 view([-171.539 51.652]) %poled dolu menší
                title("time t: "+t(i))
                pause(t(i+1)-t(i))
            end
        end
        function plotNodes(obj)
            hold off
            plot3(obj.nodes(1,:),obj.nodes(2,:),obj.nodes(3,:), '*')
            hold on
            for i = 1:obj.bars.count
                mid_point = [obj.s(6*(i-1)+1);obj.s(6*(i-1)+2);obj.s(6*(i-1)+3)];
                transform_matrix = obj.rMatrix()*obj.Tpx(obj.s(6*(i-1)+4))*obj.Tpy(obj.s(6*(i-1)+5));
                plot3(mid_point(1), mid_point(2), mid_point(3), 'o')
                point1 = mid_point+transform_matrix*[0;0;-obj.bars.lengths(i)/2;1];
                point2 = mid_point+transform_matrix*[0;0; obj.bars.lengths(i)/2;1];
                plot3([point1(1) point2(1)],[point1(2) point2(2)],[point1(3) point2(3)], 'Color','red');
            end
            transform_matrix = obj.rMatrix()*obj.Tpx(obj.s(end-2))*obj.Tpy(obj.s(end-1))*obj.Tpz(obj.s(end));
            mid_point = obj.s(end-5:end-3);
            points = zeros(3,3);
            for i = 1:3
                r_end_efector = obj.frames.radius_top*[obj.angle2vector(120*(i)+obj.frames.rotation2),0,1]';
                points(:,i) = mid_point+transform_matrix*r_end_efector;
            end
            plot3([points(1,1) points(1,2)], [points(2,1) points(2,2)], [points(3,1) points(3,2)])
            plot3([points(1,1) points(1,3)], [points(2,1) points(2,3)], [points(3,1) points(3,3)])
            plot3([points(1,3) points(1,2)], [points(2,3) points(2,2)], [points(3,3) points(3,2)])

        end
    end
    %High tear
    methods(Access = private)
        function YD = stepFK(obj, t, Y)
            obj.s = Y(1:42,1);
            obj.sd = Y(43:end,1);
            obj.stateToNodes()
            obj.matrixFK()
            obj.vectorFK(t)

            res = obj.megaMatrixFK\obj.megaRightSideFK;
            YD = [obj.sd; res(1:42)];
            %disp("Reakce: "+norm(res(43:end)))
            clc
        end
        function matrixFK(obj)
            if isempty(obj.M)
                obj.massMatrix()
            end
            obj.jacobiMatrix()
            obj.megaMatrixFK = [obj.M, obj.Phi';obj.Phi, zeros(18,18)];
        end
        function vectorFK(obj,t)
            if isempty(obj.W)
                obj.constantVector()
            end
            obj.cablesTransformationMatrix()
            obj.deritiveOfJacobiMatrix()
            obj.stringForces()
            obj.constrainResiduum()
%             disp("t: "+t)
            c1=0;
            c2=1;
            if norm(obj.residuum)<0.1
                c1 = min(1, 0.1*t^2);
            else
                disp("Oh no t: "+t)
            end
            obj.megaRightSideFK = [(c2*obj.W+c1*obj.Tc*obj.cable_forces); -(obj.PhiD*obj.sd+2*obj.alf*obj.Phi*obj.sd+obj.bet^2*obj.residuum)];
        end
        function jacobiMatrix(obj)
            %Vazbové rovnice se píší pro dva typy připojení -
            %1) rám - tyč
            %   -tj. rovnice 1-9
            %2) endEfektor - tyč
            %   -tj rovnice 10-18
            %   s = [x1 y1 z1 px1 py1 pz1 x2 ...]
            %   q = [x7 y7 z7 px7 py7 pz7]
            % Prerekvizity: nic
            % Rovnice pro rám %
            phi = zeros(3*6, 6*7);
            for i = 1:3 %i - aktuálně řešená tyč
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                inputs = {obj.s(current_vars_indexes(4)),obj.s(current_vars_indexes(5)), 0, 1, 1, 1};
                phi(current_eq_indexes, current_vars_indexes) = obj.barDeritive(inputs{:},-1,i);
            end
            %Rovnice pro endefektor
            for i = 4:6
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                %Pro tyče zůstávají
                inputs = {obj.s(current_vars_indexes(4)),obj.s(current_vars_indexes(5)), 0, 1, 1, 1};
                phi(current_eq_indexes, current_vars_indexes) = obj.barDeritive(inputs{:},1,i);
                %Pro endefektor se musí přidat extra výrazy
                phi(current_eq_indexes, end-5:end) = -phiEndEfectorDeritive(obj, (i-3));
            end
            obj.Phi = phi;
        end
        function deritiveOfJacobiMatrix(obj)
            %Derivace jakobiánu vazbových rovnice se píší pro dva typy připojení -
            %1) rám - tyč
            %   -tj. rovnice 1-9
            %2) endEfektor - tyč
            %   -tj rovnice 10-18
            %   s = [x1 y1 z1 px1 py1 pz1 x2 ...]
            %   q = [x7 y7 z7 px7 py7 pz7]
            % Prerekvizity: nic
            % Rovnice pro rám %
            phid = zeros(3*6, 6*7);
            for i = 1:3 %i - aktuálně řešená tyč
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                inputs = {obj.bars.alpha(i),...
                    obj.bars.beta(i),...
                    0,...
                    obj.sd(6*(i-1)+1),...
                    obj.sd(6*(i-1)+2),...
                    obj.sd(6*(i-1)+3),...
                    obj.sd(6*(i-1)+4),...
                    obj.sd(6*(i-1)+5),...
                    obj.sd(6*(i-1)+6)};
                phid(current_eq_indexes, current_vars_indexes) = obj.barDeritiveD(inputs{:},-1, i);
            end
            %Pro end efektor
            for i = 4:6
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                %Pro tyče zůstávají
                inputs = {obj.bars.alpha(i),...
                    obj.bars.beta(i),...
                    0,...
                    obj.sd(6*(i-1)+1),...
                    obj.sd(6*(i-1)+2),...
                    obj.sd(6*(i-1)+3),...
                    obj.sd(6*(i-1)+4),...
                    obj.sd(6*(i-1)+5),...
                    obj.sd(6*(i-1)+6)};
                phid(current_eq_indexes, current_vars_indexes) = obj.barDeritiveD(inputs{:},1,i);
                %Pro endefektor se musí přidat extra výrazy
                phid(current_eq_indexes, end-5:end) = -phiEndEfectorDeritiveD(obj, (i-3));
            end
            obj.PhiD = phid;
        end
        function cablesTransformationMatrix(obj)
            %Transformační matice pro výpočet sil na tělesa od lan
            %Prerekvizity: statesToNodes
            %Cílem je vytvořit matici Tc takovou, že platí:
            %Tc*F = G
            %Kde F = sloupcový vektor obsahující síly v lanech
            %a G je vektor silových účinků na jednotlivé tyče
            obj.Tc = zeros(42, obj.cables.count);

            %Potřebuji projet všechny tyče
            for current_bar_index = 1:obj.bars.count
                %Dolní / horní číslo uzlu:
                nodes_cols = obj.bars.from_to(current_bar_index,:);
                %Příslušné stavové veličiny
                current_s = obj.s(6*(current_bar_index-1)+(1:6));
                r = [0;0;obj.bars.lengths(current_bar_index)/2;1];
                r1 = obj.rMatrix()*obj.Tpx(-current_s(4))*obj.Tpy(current_s(5))*(-r);
                r2 = obj.rMatrix()*obj.Tpx(-current_s(4))*obj.Tpy(current_s(5))*(r);
                %Řeší se zde spodní část tyče - najdi všechny lana, co se
                %dotýkají spodního uzlu tj nodes_cols(1)
                [cable_numbers, cable_from_to_index] = find(obj.cables.from_to == nodes_cols(1));
                for current_cable_index = 1:numel(cable_numbers)
                    swap_index = [2,1];
                    %z uzlu n1 jde lano do uzlu n2
                    cable_number = cable_numbers(current_cable_index);
                    n1 = obj.nodes(:,nodes_cols(1));
                    n2 = obj.nodes(:,obj.cables.from_to(cable_number,swap_index(cable_from_to_index(current_cable_index))));

                    %Vektor mezi uzly
                    obj.Tc(6*(current_bar_index-1)+(1:3),cable_number) = (n1-n2)/norm(n1-n2);
                    obj.Tc(6*(current_bar_index-1)+(4:6),cable_number) = cross(r1,obj.Tc(6*(current_bar_index-1)+(1:3),cable_number));
                end
                %Řeší se zde horní část tyče
                [cable_numbers, cable_from_to_index] = find(obj.cables.from_to == nodes_cols(2));
                for current_cable_index = 1:numel(cable_numbers)
                    swap_index = [2,1];
                    %z uzlu n1 jde lano do uzlu n2
                    cable_number = cable_numbers(current_cable_index);
                    n1 = obj.nodes(:,nodes_cols(2));
                    n2 = obj.nodes(:,obj.cables.from_to(cable_number,swap_index(cable_from_to_index(current_cable_index))));

                    obj.Tc(6*(current_bar_index-1)+(1:3),cable_number) = (n1-n2)/norm(n1-n2);
                    obj.Tc(6*(current_bar_index-1)+(4:6),cable_number) = cross(r2, obj.Tc(6*(current_bar_index-1)+(1:3),cable_number));
                end
            end
        end
        function stringForces(obj)
            %Vypočítá absolutní velikost síly v lanech
            %Prerekvizity: stateToNodes, calculateNodesVelocities
            %Z toho co je na generátoru vychází, že délky paralelních
            %pružin jsou větší, resp normální volné délky jsou o 5% kratší
            obj.nodesVelocities()
            l = sqrt(sum((obj.nodes*obj.cables.connectivity_matrix').^2));
            obj.cable_forces = zeros(18,1);
            for i = 1 : 18
                li = l(i);
                l0i = obj.cables.free_length(i);
                l0pi = obj.cables.free_length(i);
                dli = (li-l0i);
                dlpi = (li-l0pi);
                ksi = obj.cables.specific_stiffness(i);
                Kpi = obj.cables.stiffness(i);
                vi = 0;
                bsi = obj.cables.specific_dumpings(i);
                Bpi = obj.cables.dumpings(i);
                if dli < 0; dli = 0; end
                if dlpi < 0; dlpi = 0; end
                obj.cable_forces(i) = -(ksi*dli/l0i+Kpi*dlpi+bsi/l0i*vi+Bpi*vi);
            end
        end
        function constrainResiduum(obj)
            residuum = zeros(3*6,1);
            for i = 1:3 %i - aktuálně řešená tyč
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                inputs = obj.s(current_vars_indexes);
                ri = [obj.frames.radius_bot*obj.angle2vector((i-1)*120), 0]';
                residuum(current_eq_indexes) =[inputs(1:3)] + obj.rMatrix()*obj.Tpx(inputs(4))*obj.Tpy(inputs(5))*obj.Tpz(inputs(6))*[0;0;-obj.bars.lengths(i)/2;1]...
                    -ri;
            end
            for i = 4:6 %i - aktuálně řešená tyč
                current_vars_indexes = 6*(i-1)+(1:6);
                current_eq_indexes = 3*(i-1)+(1:3);
                inputs = obj.s(current_vars_indexes);
                r_end_efector = obj.frames.radius_top*[obj.angle2vector(120*(i)+obj.frames.rotation2),0,1]';
                residuum(current_eq_indexes) =[inputs(1:3)] + obj.rMatrix()*obj.Tpx(inputs(4))*obj.Tpy(inputs(5))*obj.Tpz(inputs(6))*[0;0;obj.bars.lengths(i)/2;1]...
                    -(obj.s(end-5:end-3)+obj.rMatrix()*obj.Tpx(obj.s(end-2))*obj.Tpy(obj.s(end-1))*obj.Tpz(obj.s(end))*r_end_efector);
            end
            obj.residuum = residuum;
        end
        function nodesVelocities(obj)
            obj.nodes_velocity = zeros(3,obj.frames.nodes_count);
            velocities = zeros(3, obj.frames.nodes_count);
            for i = 1:6
                current_vars_indexes = 6*(i-1)+(1:6);
                inputs = {obj.sd(current_vars_indexes(1)),...
                    obj.sd(current_vars_indexes(2)),...
                    obj.sd(current_vars_indexes(3)),...
                    obj.s(current_vars_indexes(4)),obj.s(current_vars_indexes(5)), 0,...
                    obj.sd(current_vars_indexes(4)),obj.sd(current_vars_indexes(5)), 0};
                %Spodní části tyče
                velocities(:, obj.bars.from_to(i,1)) = obj.nodeVelocity(inputs{:},-1, i);
                %Horní části tyče
                velocities(:, obj.bars.from_to(i,2)) = obj.nodeVelocity(inputs{:}, 1, i);
            end
        end
    end
    %Mid tear
    methods(Access = private)
        function phi = barDeritive(obj, phix, phiy, phiz, phixd, phiyd, phizd, dir, bar_number)
            %Vrátí gradient polohy koncového bodu tyče - grad(midPoint+TpxTpyTpz*r)
            r = [0;0;dir*obj.bars.lengths(bar_number)/2;1];
            phi = [eye(3),...
                obj.rMatrix()*obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.Tpz(phiz)*r,...
                obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*r,...
                obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*r];
        end
        function phi = phiEndEfectorDeritive(obj, point_number)
            %Point number - číslo rohu range 1-3
            phix = obj.s(end-2);
            phiy = obj.s(end-1);
            phiz = obj.s(end);
            r = obj.frames.radius_top*[obj.angle2vector(120*(point_number)+obj.frames.rotation2),0,1]';
            phi = [eye(3),...
                obj.rMatrix()*obj.Tpx(phix)*obj.DTpx(1)*obj.Tpy(phiy)*obj.Tpz(phiz)*r,...
                obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(1)*obj.Tpz(phiz)*r,...
                obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(1)*r];
        end
        function phiD = barDeritiveD(obj, phix, phiy, phiz, xd, yd, zd, phixd, phiyd, phizd, dir, bar_number)
            %phizd neimplementováno
            r = [0;0;dir*obj.bars.lengths(bar_number)/2;1];
            phiD4 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.DTpx(1)*obj.Tpy(phiy)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.DTpx(1)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.DTpx(1)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            phiD5 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.DTpy(1)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.DTpy(1)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(1)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            phiD6 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(1)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*obj.DTpz(1)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*obj.DTpz(1)*r;
            phiD = [zeros(3), phiD4(1:3),phiD5(1:3), phiD6(1:3)];
        end
        function phiD = phiEndEfectorDeritiveD(obj, point_number)
            %Point number - číslo rohu range 1-3
            r = obj.frames.radius_top*[obj.angle2vector(120*(point_number)+obj.frames.rotation2),0,1]';
            phix = obj.s(end-2);
            phiy = obj.s(end-1);
            phiz = obj.s(end);
            phixd = obj.sd(end-2);
            phiyd = obj.sd(end-1);
            phizd = obj.sd(end);
            phiD4 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.DTpx(1)*obj.Tpy(phiy)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.DTpx(1)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.DTpx(1)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            phiD5 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.DTpy(1)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.DTpy(1)*obj.Tpz(phiz)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(1)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            phiD6 = obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(1)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*obj.DTpz(1)*r+...
                obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*obj.DTpz(1)*r;
            phiD = [zeros(3), phiD4(1:3),phiD5(1:3), phiD6(1:3)];
            obj.PhiD = phiD;
        end
        function velocity = nodeVelocity(obj, vx, vy, vz, phix, phiy, phiz, phixd, phiyd, phizd, dir, bar_number)
            %Vrátí rychlost uzlu
            r = [0;0;dir*obj.bars.lengths(bar_number)/2;1];
            velocity_from_x_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.DTpx(phixd)*obj.Tpy(phiy)*obj.Tpz(phiz)*r;
            velocity_from_y_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*r;
            velocity_from_z_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            velocity = [vx;vy;vz]+velocity_from_x_rotation+velocity_from_y_rotation+velocity_from_z_rotation;
        end

        %Konstantní matice / vektory
        function massMatrix(obj)
            %Vypočítá matici hmotnosti
            %Prerekvizity: nic
            vector = zeros(1,(obj.bars.count+1)*6);
            for i = 1 : obj.bars.count
                vector(6*(i-1)+1:(6*(i-1))+6) = [obj.bars.masses(i)*ones(1,3), obj.bars.inertias_x(i), obj.bars.inertias_y(i), obj.bars.inertias_z(i)];
                obj.M = diag(vector);
            end
            obj.M((end-5):end,(end-5):end) = diag([obj.bars.masses(i)*3*ones(3,1); 0.01*ones(3,1)]);
        end
        function constantVector(obj)
            %Vypočítá tíhové síly
            %Pro horní platformu je hmotnost dána trojnásobkem tyčí
            %Prerekvizity: nic
            obj.W = zeros(6,7);
            obj.W(3,:) = [reshape(obj.bars.masses,1,[]), obj.bars.masses(1)*3]*(obj.g);
            obj.W = obj.W(:);
        end

        %Přepočty
        function initialConditionsToStates(obj)
            obj.s = zeros(6*(obj.bars.count+1),1);
            obj.sd = zeros(6*(obj.bars.count+1),1);
            %Počáteční podmínky pro tyče
            for i = 1:obj.bars.count
                obj.s(6*(i-1)+(1:6),1) = [reshape(obj.bars.mid_points(:,i),[],1);
                    obj.bars.alpha(i);
                    obj.bars.beta(i);
                    0];
                obj.sd((i-1)+(1:6),1) = zeros(6,1);
            end
            i = 7;
            %Počáteční podmínky pro end efektor
            obj.s(6*(i-1)+(1:6),1) = [obj.frames.x;obj.frames.y;obj.frames.z;
                obj.frames.px;obj.frames.py;obj.frames.pz];
            obj.sd(6*(i-1)+(1:6),1) = [obj.frames.xd;obj.frames.yd;obj.frames.zd;
                obj.frames.pxd;obj.frames.pyd;obj.frames.pzd];
        end
        function stateToNodes(obj)
            for current_bar_index = 1:obj.bars.count
                nodes_cols = obj.bars.from_to(current_bar_index,:);
                current_s = obj.s(6*(current_bar_index-1)+(1:6));
                position = current_s(1:3);
                r = [0;0;obj.bars.lengths(current_bar_index)/2;1];
                obj.nodes(:,nodes_cols(1)) = position+obj.rMatrix()*obj.Tpx(current_s(4))*obj.Tpy(current_s(5))*(-r);
                obj.nodes(:,nodes_cols(2)) = position+obj.rMatrix()*obj.Tpx(current_s(4))*obj.Tpy(current_s(5))*(r);
            end
        end
    end
    %Low tear
    methods(Access=private)
        %Transformační matice
        function Tpx = Tpx(~, phi_y)
            Tpx = [ 1 0 0 0;...
                0 cos(phi_y) -sin(phi_y) 0;...
                0 sin(phi_y) cos(phi_y) 0;...
                0 0 0 1];
        end
        function Tpy = Tpy(~, phi_y)
            Tpy = [ cos(phi_y) 0 sin(phi_y) 0;...
                0 1 0 0;...
                -sin(phi_y) 0 cos(phi_y) 0;...
                0 0 0 1];
        end
        function Tpz = Tpz(~, phi_z)
            Tpz = [ cos(phi_z) -sin(phi_z) 0 0;...
                sin(phi_z) cos(phi_z) 0 0;...
                0 0 1 0;...
                0 0 0 1];
        end
        function Tx = Tx(~, x)
            Tx = eye(4);
            Tx(1,4) = x;
        end
        function Ty = Ty(~, y)
            Ty = eye(4);
            Ty(2,4) = y;
        end
        function Tz = Tz(~, z)
            Tz = eye(4);
            Tz(3,4) = z;
        end

        %Matice pro derivování
        function DTpx = DTpx(~, dx)
            DTpx = zeros(4);
            DTpx(2,3) = -dx;
            DTpx(3,2) = dx;
        end
        function DTpy = DTpy(~, dy)
            DTpy = zeros(4);
            DTpy(1,3) = dy;
            DTpy(3,1) = -dy;
        end
        function DTpz = DTpz(~, dz)
            DTpz = zeros(4);
            DTpz(1,2) = -dz;
            DTpz(2,1) = dz;
        end

        %Ostatní
        function rM = rMatrix(~)
            %Vrátí matici, která redukuje rozměr
            %např rozšířené vektory r jsou 4x1, po aplikaci
            %rM*r získám vektor 3x1 - ořezaný o 1 na konci
            rM = [eye(3), zeros(3,1)];
        end
    end
end

