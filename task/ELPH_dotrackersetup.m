%%% Wechselt in den Kalibrationsmodus und führt die Kalibration durch

function ELPH_dotrackersetup(el)

Eyelink( 'Command', 'heuristic_filter = ON');
ListenChar(2);
Eyelink( 'StartSetup',1);		% start setup mode
ListenChar(1);
Eyelink( 'WaitForModeReady', el.waitformodereadytime );

return	
