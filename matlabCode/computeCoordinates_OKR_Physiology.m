%% computeCoordinates_OKR_Physiology
% Figures out the polar coordinates of each cell and saves it to the data
% structure

cells = fieldnames(data);

radius = data.centerData.transformedCoordinates.medianRadius; %in um
try
    northPoint = data.centerData.transformedCoordinates.northPoint;
catch
    northPoint = [0 0];
end
[northPoint_angle, ~] = cart2pol(northPoint(1), northPoint(2));
angleOffset = northPoint_angle - pi/2; %this is the offset of the tissue. Always add it to the theta coordinate of a cell

if isfield(data, 'centerData')
    for i = 1:size(cells, 1)
        cell_i = cells{i, 1};
        if ~strcmp(cell_i, 'centerData') %one field is centerData, don't want to look at this one here

            %check whether there is a metadata field called coordinates for
            %each cell. Case insensitive.
            metaFields = fieldnames(data.(cell_i).epochs(1).meta);
            if ismember('Coordinates', metaFields)
                cord_fieldname = 'Coordinates';
            elseif ismember('coordinates', metaFields')
                cord_fieldname = 'coordinates';
            else
                cord_fieldname = '';
            end

            added = 0;
            if cord_fieldname %if there's coordinate data for this cell
                cartesian_cords = data.(cell_i).epochs(1).meta.(cord_fieldname);

                %strip spaces
                cartesian_cords = cartesian_cords(~isspace(cartesian_cords));

                %turn string into a cell array
                cartesian_cords_array = split(cartesian_cords, ',');

                if numel(cartesian_cords_array) == 2
                    added = 1;
                    xCord = str2num(cartesian_cords_array{1});
                    yCord = str2num(cartesian_cords_array{2});


                    [t, r] = cart2pol(xCord, yCord);
                    t = angleOffset + t;

                    r_normalized = r/radius;

                    data.(cell_i).coordinates.polar = [t, r_normalized];
                    data.(cell_i).coordinates.polar_info = 'Polar coordinates in terms of [theta, rho] where rho is normalized to the radius of the retina';
                    data.(cell_i).coordinates.cartesian_uncorrected = [xCord, yCord];
                    data.(cell_i).coordinates.centerData = data.centerData;
                end
            end

            if ~added
                if isfield(data.(cell_i), 'coordinates')
                    warning(['Did not add coordinates or center data for, ', cell_i, 'but previous data found'])
                else
                    data.(cell_i).coordinates = 'unknown';
                    warning(['Coordinates unknown for', cell_i, 'added as unknown'])
                end
            end
        end
    end
    
else %if centerData is not included
    warning('Center data not included')
end

clear added angleOffset cartesian_cords cartesian_cords_array cell_i cells...
    cord_fieldname i metaFields northPoint northPoint_angle northPoint_polar...
    northPoint_polar_distance r r_normalized radius t t_normalized xCord yCord