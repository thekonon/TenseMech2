classdef TensegritySettings
    %Třída obsahuje základní parametry úlohy
    properties (Access = protected)
        bars    bars
        cables  cables
        frames  frames
        nodes
    end

    methods
        function obj = TensegritySettings(varargin)
            %% Nastaví hodnoty společné pro všechny úlohy
            %Jedná se o vyplnění objektů tyčí, lan a rámů
%             setProperties(obj,nargin,varargin{:})
            obj.inicialization();
            %% Možnost přepsání hodnot, které jsou například optimalizovány
            %Takto by vypadalo přepsání tuhosti lana
            %první argument - index lana
            %druhý - jeho měrná tuhost
            %obj.cables.specific_stiffness(varargin{1}) = varargin{2};
            %% Výpočet matic konektivity
            obj.calculateConnectivityMatrixes()
            %% Zde se ze zadaných hodnot dopočítají polohy uzlů
            obj.initialNodesGuess();
            
            %Délky lan je třeba dopočítat z původního odhadu (délky
            %paralelních pružin například)
            obj.cables.free_length = sqrt(sum((obj.nodes*obj.cables.connectivity_matrix').^2));
%             load steady_state_nodes.mat nodes
%             obj.nodes = nodes'; %Které je možné přepsat třeba přesnějším
%                                     řešením
            %% Dopočet počátečních souřadnic
            %tj natočení tyčí, polohy středů pro SimScape
            obj.calculateInitialCoordinates();
            %% Dopočet natočení end-efektoru
            obj.calculateEndEfector();

        end
    end
    methods(Access = private)
        function inicialization(obj)
            %Nastaví pevné parametry společné pro všechny úlohy
            obj.setBars();
            obj.setCables();
            obj.setFrames();
        end
        %Inicialization
        function setBars(obj)
            %vlastnosti tyčí jsou uloženy v instanci objektu bars
            obj.bars = bars();
            obj.bars.count = 6;     %Počet tyčí
            obj.bars.lengths =      ones(obj.bars.count,1)*0.3;         %Délky
            obj.bars.diameters =    ones(obj.bars.count,1)*0.01;        %Průměry
            obj.bars.densitys =     ones(obj.bars.count,1)*3000;        %Hustoty
            %Hmotnové vlastnosti
            obj.bars.masses =   obj.bars.lengths.*obj.bars.diameters.^2/4.*obj.bars.densitys;
            obj.bars.inertias_x = 1/12*obj.bars.masses.*obj.bars.lengths.^2;
            obj.bars.inertias_y = obj.bars.inertias_x;
            obj.bars.inertias_z = obj.bars.masses.*(obj.bars.diameters/2).^2;
            %Tabulka z do - pro vytvoření tabulky konektivity
            obj.bars.from_to = [1 4;...
                2 5
                3 6
                7 11
                8 12
                9 10];
            %Rozšíření - z historických důvodů, pokud se modeluje horní
            %platforma jako tři tyče
            obj.bars.from_to_extended = [1 4;...
                2 5
                3 6
                7 11
                8 12
                9 10
                10 11
                11 12
                10 12];
        end
        function setCables(obj)
            obj.cables = cables();
            obj.cables.count = 18;  %Počet lan
            obj.cables.specific_stiffness = ones(obj.cables.count,1)*200;
            obj.cables.specific_dumpings = ones(obj.cables.count,1)*30;
            obj.cables.stiffness = ones(obj.cables.count,1)*25;
            obj.cables.dumpings = ones(obj.cables.count,1)*10;
            %varialbe cables - zde se nastavuje, která lana jsou řízena
            obj.cables.variable_cables_indexes = [7,8,9,13,14,15];
            obj.cables.fixed_cables_indexes = setdiff(1:obj.cables.count, obj.cables.variable_cables_indexes,'sorted');
            %Opět tabulka z do
            obj.cables.from_to = ...
                [1 5    %1
                1 7     %2
                2 6     %3
                2 8     %4
                3 4     %5
                3 9     %6
                4 7    %7
                4 9     %8
                5 7     %9
                5 8     %10
                6 8     %11
                6 9    %12
                6 12    %13
                9 12     %14
                4 10     %15
                7 10    %16
                5 11    %17
                8 11];  %18
        end
        function setFrames(obj)
            %Vlastnosti rámů
            obj.frames = frames();
            obj.frames.radius_mid = 0.12;       %velikost kružnice středu
            obj.frames.rotation1 = -120;        %O kolik je pootočena prostřední platforma vůči spodku
            obj.frames.rotation2 = -60;         %O kolik je pootočena horní platforma
            %Počet uzlů:
            obj.frames.nodes_count = max([obj.bars.from_to(:);obj.cables.from_to(:)]);
            %Pevné uzly jsou zde hlavně z historických důvodů
            obj.frames.fixed_nodes = zeros(1,obj.frames.nodes_count);
            obj.frames.fixed_nodes([1,2,3]) = [1,1,1];
            %Layers shift je posun mezi dolní a horní částí tensegrity
            obj.frames.layers_shift = 0.05;
            %Zde jsou jen tloušťky platforem
            obj.frames.width_bot = 0.005;
            obj.frames.width_top = 0.005;
            %A poloměry základen
            obj.frames.radius_bot = 0.1;
            obj.frames.radius_top = 0.1;
            obj.frames.type = 1;
        end
        function initialNodesGuess(obj)
            %Do obj.nodes vyplňí počáteční polohu uzlů
            obj.nodes = zeros(3,max(obj.bars.from_to(:)));
            %Bot nodes
            obj.nodes(1:2,1) = obj.frames.radius_bot*obj.angle2vector(0);
            obj.nodes(1:2,2) = obj.frames.radius_bot*obj.angle2vector(120);
            obj.nodes(1:2,3) = obj.frames.radius_bot*obj.angle2vector(240);
            %Mid nodes
            obj.nodes(1:2,4) = obj.frames.radius_mid*obj.angle2vector(0+obj.frames.rotation1);
            obj.nodes(1:2,5) = obj.frames.radius_mid*obj.angle2vector(120+obj.frames.rotation1);
            obj.nodes(1:2,6) = obj.frames.radius_mid*obj.angle2vector(240+obj.frames.rotation1);
            obj.nodes(1:2,7) = obj.frames.radius_mid*obj.angle2vector(0+obj.frames.rotation1+60);
            obj.nodes(1:2,8) = obj.frames.radius_mid*obj.angle2vector(120+obj.frames.rotation1+60);
            obj.nodes(1:2,9) = obj.frames.radius_mid*obj.angle2vector(240+obj.frames.rotation1+60);
            %Top nodes
            obj.nodes(1:2,10) = obj.frames.radius_top*obj.angle2vector(0+obj.frames.rotation2);
            obj.nodes(1:2,11) = obj.frames.radius_top*obj.angle2vector(120+obj.frames.rotation2);
            obj.nodes(1:2,12) = obj.frames.radius_top*obj.angle2vector(240+obj.frames.rotation2);

            height_first_floor = sqrt(obj.bars.lengths(1)^2 - sum(diff(obj.nodes(1:2,obj.bars.from_to(1,:))').^2));
            height_second_floor = sqrt(obj.bars.lengths(1)^2 - sum(diff(obj.nodes(1:2,obj.bars.from_to(end,:))').^2));

            obj.nodes(3,:) = repelem([0,...     první patro
                height_first_floor,...      druhé patro
                height_first_floor-obj.frames.layers_shift,...  podpatro
                height_first_floor+height_second_floor-obj.frames.layers_shift], 3);    %Vrchní patro
        end
        function calculateInitialCoordinates(obj)
            %Středy tyčí
            coresponding_nodes_bars = obj.nodes(:,obj.bars.from_to');
            obj.bars.mid_points = (coresponding_nodes_bars(:,1:2:end)+coresponding_nodes_bars(:,2:2:end))/2;

            %Natočení tyčí
            obj.bars.rotation_matrixes = zeros(3,3,obj.bars.count);
            obj.bars.alpha = zeros(obj.bars.count,1);
            obj.bars.beta = zeros(obj.bars.count,1);
            for i = 1:obj.bars.count
                [alpha,beta,rotation_matrix] = obj.transformationMatrix(coresponding_nodes_bars(:, 2*i-1),coresponding_nodes_bars(:, 2*i));
                obj.bars.alpha(i) = alpha;
                obj.bars.beta(i) = beta;
                obj.bars.rotation_matrixes(:,:,i) = rotation_matrix;
            end
        end
        function calculateEndEfector(obj)
            temp_nodes = obj.nodes(:,end-2:end);
            obj.frames.mid_point = mean(temp_nodes,2);
            obj.frames.x = obj.frames.mid_point(1);
            obj.frames.y = obj.frames.mid_point(2);
            obj.frames.z = obj.frames.mid_point(3);
            obj.frames.xd = 0;
            obj.frames.yd = 0;
            obj.frames.zd = 0;


            %abc koeficienty rovnice ax+by+cz+1 = 0
            abc = (temp_nodes')\(ones(3,1));
            obj.frames.px = atan2(abc(1), abc(3));
            obj.frames.py = atan2(abc(2), abc(3));
            obj.frames.pz = 0; 
            obj.frames.pxd = 0;
            obj.frames.pyd = 0;
            obj.frames.pzd = 0;
        end
    end
    methods(Access = protected) %Pomocné funkce
        function vector = angle2vector(~, angle_d)
            vector = [cosd(angle_d), sind(angle_d)];
        end
        function [alpha, beta, R] = transformationMatrix(~, N1, N2)
            dr = N2-N1;
            dx = dr(1);
            dy = dr(2);
            dz = dr(3);
            Rx = @(phi)[1 0 0; 0 cos(phi) -sin(phi); 0 sin(phi) cos(phi)];
            Ry = @(phi)[cos(phi) 0 sin(phi); 0 1 0; -sin(phi) 0 cos(phi)];

            alpha = -atan2(dy,dz);
            beta = atan2(dx,sqrt(dy.^2+dz.^2));

            R1 = Rx(alpha);
            R2 = Ry(beta);
            R = R1*R2;
        end

        function calculateConnectivityMatrixes(obj)
            %Matice konektivity pro tyče
            n = size(obj.bars.from_to,1);       %Celkem prvků;
            m = obj.frames.nodes_count;         %Celkem uzlů;
            C = zeros(n,m);
            for i = 1:n
                C(i,:) = obj.e(obj.bars.from_to(i,2), m) - obj.e(obj.bars.from_to(i,1), m);
            end
            obj.bars.connectivity_matrix = C;
            obj.bars.mid_points_connectivity_matrix = abs(C);
            %Matice konektivity pro lana
            n = size(obj.cables.from_to,1);  %Celkem prvků;
            m = obj.frames.nodes_count;      %Celkem uzlů;
            C = zeros(n,m);
            for i = 1:n
                C(i,:) = obj.e(obj.cables.from_to(i,2), m) - obj.e(obj.cables.from_to(i,1), m);
            end
            obj.cables.connectivity_matrix = C;
        end
        function vec = e(~,i,n)
            vec = zeros(n,1);
            vec(i) = 1;
        end
    end
end