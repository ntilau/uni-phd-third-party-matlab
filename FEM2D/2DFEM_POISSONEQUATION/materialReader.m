function material = materialReader( dateiName ) 

% erzeugung der Dateinamen 
matName = strcat( dateiName, '.mat' );

% �ffene der Datei
fid = fopen( matName );
if fid == -1
    error( 'kann %s nicht �ffnen', matName );
end
matDim = fscanf( fid, '%i' , 1 );

% lesen der Materialien
for cnt = 1:matDim
    matId = fscanf( fid, '%i' , 1 );
    if matId ~= cnt
        error( 'falschese Materialformat' );
    end
    str = fscanf( fid, '%s' , 1 );
    if strcmp(str,'epsilon') == false
        error( 'falschese Materialformat' );
    end
    epsRel = fscanf( fid, '%f' , 1 );
    str = fscanf( fid, '%s' , 1 );
    if strcmp(str,'mu') == false
        error( 'falschese Materialformat' );
    end
    muRel = fscanf( fid, '%f' , 1 );
    material( matId ).epsilonRelativ = epsRel;
    material( matId  ).muRelativ = muRel;
end
fclose( fid );

