def parse_storm_file(filename):
    storm_list = []
    current_storm = None
    
    with open(filename, 'r') as file:
        for line in file:
            # Remove leading/trailing whitespaces
            line = line.strip()
            
            if line.startswith('start'):
                # If we encounter 'start', we finish the current storm dict and start a new one
                if current_storm is not None:
                    storm_list.append(current_storm)
                
                # Initialize a new storm dictionary
                parts = line.split()
                current_storm = {
                    'lat': [],
                    'lon': [],
                    'slp': [],
                    'wind': [],
                    'surface_pressure': [],
                    'timestamp': []
                }
                continue
            
            # For non-'start' lines, read the data points
            parts = line.split()
            
            # Extract the values for the storm dict
            lon = float(parts[1])
            lat = float(parts[2])
            slp = float(parts[3])
            wind = float(parts[4])
            surface_pressure = float(parts[5])
            year = int(parts[6])
            month = int(parts[7])
            day = int(parts[8])
            hour = int(parts[9])
            
            # Append to the current storm data
            current_storm['lon'].append(lon)
            current_storm['lat'].append(lat)
            current_storm['slp'].append(slp)
            current_storm['wind'].append(wind)
            current_storm['surface_pressure'].append(surface_pressure)
            current_storm['timestamp'].append(f'{year}-{month:02d}-{day:02d} {hour:02d}:00')
        
        # Append the last storm
        if current_storm is not None:
            storm_list.append(current_storm)
    
    return storm_list
