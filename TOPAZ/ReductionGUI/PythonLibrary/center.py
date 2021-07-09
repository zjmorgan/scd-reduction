def center(h, k, l, center_type):
    """ Function to test for allowed (True) and not allowed (False)
    peaks due to centering."""
    
    if center_type == 'P':
        return True
    
    if center_type == 'A':
        sum = k + l
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'B':
        sum = h + l
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'C':
        sum = h + k
        if (sum % 2) == 0: return True
        return False
        
    if center_type == 'F':
        sum = h + k
        if (sum % 2) != 0: return False
        sum = h + l
        if (sum % 2) != 0: return False
        sum = k + l
        if (sum % 2) != 0: return False
        return True
        
    if center_type == 'I':
        sum = h + k + l
        if (sum % 2) != 0: return False
        return True
        
    if center_type == 'R':
        sum = -h + k + l
        if (sum % 3) != 0: return False
        return True
        
    print('Centering type not P, A, B, C, F, I or R.')
    
# test
if __name__ == '__main__':
    h = 1
    k = 5
    l = 3
    print(h,k,l)

    center_type = 'F'
    print(center_type)

    x = center(h, k, l, center_type)

    if x == True:
        print('x is True')
    if x == False:
        print('x is False')

        