import numpy as np 

def df2info(sc_name,center,df,max_control_points):
    """
    Extract information from a specific dataframe and generate a dictionary.
    """
    num_stars = len(df)
    if max_control_points is None or max_control_points > num_stars: 
        max_control_points = num_stars
    else:
        df = df[:max_control_points]

    ra,dec = df['ra'].values,df['dec'].values
    radec = np.stack([ra,dec]).T
    mag = df['mag'].values
    info = {'sc_name':sc_name,'num_stars':num_stars,'df':df,'center':center,'radec':radec,'mag':mag,'max_control_points':max_control_points}   

    return info