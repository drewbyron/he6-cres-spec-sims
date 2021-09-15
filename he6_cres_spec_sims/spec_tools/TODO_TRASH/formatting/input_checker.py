def isFloat(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def promptUser(prompt, rebuke = "", initialPrompt = True):
        
        if rebuke == "":
            rebuke = "Please type 'y' or 'n'"
        
        if (initialPrompt == True):
            print()
            choice = input(prompt)

        while (choice !="y") and (choice !="n"):
            print(rebuke)
            print()
            choice = input(prompt)
        
        return choice
        
def promptRange(prompt, rebuke, test="NONE", min=0, max=0):
        
    if (test == "LoG"):
        choice = input(prompt)
        
        while (not isFloat(choice)) or (float(choice) < min) or (float(choice) > max):
            print(rebuke)
            print()
            choice = input(prompt)
            
    if (test == "GoE"):
           choice = input(prompt)
           
           while (not isFloat(choice)) or (float(choice) < min):
               print(rebuke)
               print()
               choice = input(prompt)
               
    if (test == "LoE"):
           choice = input(prompt)
           
           while (not isFloat(choice)) or (float(choice) > max):
               print(rebuke)
               print()
               choice = input(prompt)
               
    if (test == "NONE"):
           choice = input(prompt)
           
           while (not isFloat(choice)):
               print(rebuke)
               print()
               choice = input(prompt)
               
    return float(choice)
        
        
    
