classdef tenseMech<TensegritySettings    
    %Veřejné metody
    properties(Access = public)
        Property1
    end
    %Privátní metody
    properties(Access = private)
        M
        Phi
        PhiD
        W

        %Stavové proměnné
        s
        sd
    end
    %Konstanty
    properties(Access = public,Constant)
        time_stop = 0.1
        alf = 1
        bet = 1
        g = -9.81
    end
    
    methods(Access = public)
        function obj = tenseMech()
            obj@TensegritySettings()
            old_nodes = obj.nodes;
            obj.initialConditionsToStates()
            obj.stateToNodes()
        end
    end
    %High tear
    methods(Access = private)
        
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
            phiD = [diag([xd,yd,zd]), phiD4(1:3),phiD5(1:3), phiD6(1:3)];
        end
        function phiD = phiEndEfectorDeritiveD(obj, point_number)
            %Point number - číslo rohu range 1-3
            r = obj.frames.radius_top*[obj.angle2vector(120*(point_number)+obj.frames.rotation2),0,1]';
            phix = obj.s(end-2);
            phiy = obj.s(end-1);
            phiz = obj.s(end);
            xd = obj.sd(end-5);
            yd = obj.sd(end-4);
            zd = obj.sd(end-3);
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
            phiD = [diag([xd,yd,zd]), phiD4(1:3),phiD5(1:3), phiD6(1:3)];
            obj.PhiD = phiD;
        end
        function velocity = nodeVelocity(obj, vx, vy, vz, phix, phiy, phiz, phixd, phiyd, phizd, dir, bar_number)
            %Vrátí rychlost uzlu
            r = [0;0;dir*obj.bars.lengths(bar_number)/2;1];
            phix = -phix;
            velocity_from_x_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.DTpx(-phixd)*obj.Tpy(phiy)*obj.Tpz(phiz)*r;
            velocity_from_y_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.DTpy(phiyd)*obj.Tpz(phiz)*r;
            velocity_from_z_rotation = obj.rMatrix()*obj.Tpx(phix)*obj.Tpy(phiy)*obj.Tpz(phiz)*obj.DTpz(phizd)*r;
            velocity = [vx;vy;vz]+velocity_from_x_rotation+velocity_from_y_rotation+velocity_from_z_rotation;
        end
        %přepočty
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

