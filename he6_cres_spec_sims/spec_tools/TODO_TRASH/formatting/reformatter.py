def config_tick_reformatter(dec_places=1):
    """
    Constructs a tick reformatter function to read floats as <num> x 10^{<exp>} where <num> is approximated to the given number of decimal places.
    """

    format_style = ".{}e".format(dec_places)

    def tick_reformatter(tick_val, pos):

        tick_formatted = ("{:"+format_style+"}").format(tick_val)
        tick_formatted = str(tick_formatted)

        index_of_e = tick_formatted.find('e')
        tick_power = '{0:d}'.format(int(tick_formatted[(index_of_e + 1):]))
        if tick_formatted[0] == '0':
            tick_formatted = '0'
        elif tick_power == "0":
            tick_formatted = tick_formatted[0:(index_of_e)]
        elif tick_power == "1":
        
            index_of_decimal = tick_formatted.find('.')
            if tick_formatted[index_of_decimal+2:(index_of_e)] == "":
                tick_formatted = tick_formatted[:index_of_decimal] + tick_formatted[index_of_decimal+1] + ".0"
            else:
                tick_formatted = tick_formatted[:index_of_decimal] + tick_formatted[index_of_decimal + 1] + "." + tick_formatted[index_of_decimal + 2:(index_of_e)]
        elif tick_power == "-1":
        
            if tick_formatted[0] == "-":
                tick_formatted = "-0." + tick_formatted[1] +tick_formatted[3:(index_of_e)]
            else:
                tick_formatted = "0." + tick_formatted[0] +tick_formatted[2:(index_of_e)]
            
            if tick_formatted[-1] == "0":
                tick_formatted = tick_formatted[:-1]
        else:
            tick_formatted = tick_formatted[0:(index_of_e)] + r"$\times $" +  "$10^{" + tick_power + "}$"
        return tick_formatted
    return tick_reformatter

def sci_reformatter(val,dec_places=2):
    """
    Formats floats to strings that read as <num> x 10^{<exp>} where <num> is approximated to the given number of decimal places.
    """
    format_style = ".{}e".format(dec_places)
    val_txt = ("{:" + format_style + "}").format(val)
    index_of_e = val_txt.find('e')

    val_power = '{0:d}'.format(int(val_txt[(index_of_e + 1):]))

    if val_txt[0] == '0':
        val_txt = '0'
    elif int(val_power) > -2 and int(val_power) < 0:
        format_style = ".{}f".format(abs(int(val_power)) + dec_places)
        val_txt = ("{:" + format_style + "}").format(val)
    elif int(val_power) >= 0 and int(val_power) < 2:
        format_style = ".{}f".format(dec_places)
        val_txt = ("{:" + format_style + "}").format(val)
    else:
        val_txt = val_txt[0:(index_of_e)] + r"$\times $" +  "$10^{" + val_power + "}$"
    return val_txt
    
def sigfig(val,places=2):
    """
    Calculates significant figures to number of places and returns float.
    """
    format_style = ".{}e".format(places + 1)
    val_txt = ("{:" + format_style + "}").format(val)
    index_of_e = val_txt.find('e')
    
    val_power = '{0:d}'.format(int(val_txt[(index_of_e + 1):]))
    
    val_figures = val_txt[:index_of_e-1]
    val_figures = round(float(val_figures),places - 1)
    #print(val_figures)
    num_val_figures = len(val_txt)
    
    final_figure = val_figures * pow(10,float(val_power))
    return final_figure
    

