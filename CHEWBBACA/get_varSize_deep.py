#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Friday January 25 2019

@author: rfm

######################
# Module description #
######################

To actually know the size of variables as lists, sets and dicts plus the objects
to which the contained pointers refer, it is needed to write a recursive function.
This module adresses that problem and calculates/estimates the memory footprint of objects.

##############################################################
# Notes about the size and assignment of variables in Python #
##############################################################

The 'sys.getsizeof' function from the 'sys' module gives a shallow value
for the size of a variable.

sys.getsizeof(5) # 28 bytes per integer
sys.getsizeof(5.3) # 24 bytes per float
sys.getsizeof('') # 49 bytes for an empty string
sys.getsizeof('12') # 51 bytes - 37 for empty string and 1 for each char
sys.getsizeof([]) # 64 bytes for an empty list
sys.getsizeof([1]) # 72 bytes - 64 for the empty list and 8 for each added object
sys.getsizeof(['2']) # 72 bytes - 64 for the empty list and 8 for each added object
# The list does not contain the objects, it just contains a 8 bytes pointer per object
sys.getsizeof(()) # 48 bytes for empty tuple
sys.getsizeof((1,'23')) # 64 bytes - 48 for the empty tuple and 8 for each pointer
sys.getsizeof(set()) # 224 bytes per set and size will not increase by adding objects
sys.getsizeof(dict()) # 240 per dict and size will not increase by adding objects

Python might keep memory that was used for some temporary condition, meaning
that it might reserve a lot of memory when it is not needed anymore because
the program needed it at some point but now it only needs a fraction of that
memory.

Assigning one variable to another (b = a) just results in the variables 
refering to the same value/object in memory. So, each variable in Python just refers
to a value in memory, working sort as a pointer in C. Python keeps track of the
references to an object and if the number of references drops to zero, the
garbage collector removes the object from memory.

Creating variables that refer to two different objects with equal value:
    
>>> L1 = [1,2,3]
>>> L2 = [1,2,3]

>>> L1 == L2
True

>>> L1 is L2
False

Creating variables that point to the same object:

>>> L1 = [1,2,3]
>>> L2 = L1

>>> L1 == L2
True

>>> L1 is L2
True

If we alter the object through one variable, it affects the result when we call
the other variable:

>>> L1.append(4)
>>> print(L2)
[1,2,3,4]

We can change the object to which a variable refers to as if we are just
switching a label from one object to another:
    
>>> x = 8
>>> y = x
>>> x = 100
>>> y
8

The variable/label 'x' that was pointing to the object/value 8 was reassigned 
to the value 100 but the variable 'y' still refers to the object/value 8 because
when it was assigned, the variable 'x' was still pointing to 8 and not to 100.

String interning:
Some variables can be created by assigning them to apparently different string 
objects but end up refering to the same object/value.

>>> a = "python"
>>> b = "python"
>>> a is b
True   # a and b refer to the same object!

Integer caching:
Python keeps a global list of all integers in the range [-5,256], allocating
266*24 = 6384 bytes for all those integers. The integers in this range all have
the same 'id' and working with their ids should be done with caution.

>>> a = 256
>>> b = 256
>>> a is b
True

If the integer value is outside the range [-5,256], this doe snot happen:
>>> a = 257
>>> b = 257
>>> a is b
False

Empty immutable objects:

Variables assigned to an empty tuple refer to the same object:
>>> a = ()
>>> b = ()
>>> a is b
True  # a and b both refer to the same object in memory

If the tuples are not empty, distinct objets are created:
>>> a = (1, )
>>> b = (1, )
>>> a == b
True
>>> a is b
False

Caution about Operator behavior:
L = [1,2,3]

L += [x] and L.append(x) modify the object that the variable L refers to because
the '+=' calls the __iadd__ method that modifies the objects/arguments in place.
The append method has this behavior.
But L = L + [4] will call the __add__ method that does not modify the arguments
in place and creates a new object that L will refer to.

Passing immutable objects:
def increment(n):
    n += 1

>>> a = 3
>>> increment(a)
>>> print(a)
a = 3   # a is still referring to the same object

The variable 'n' created with the function will refer to the same object as the
variable 'a' but it will be reassigned inside the function to the value 4. The
variable 'a' remains assigned to the value 3.

If we want to reassign variable 'a', we have to reassign it with the function:
def increment(n):
    n += 1
    return n

>>> a = 3
>>> a = increment(a)  # the return value of increment() is captured!
>>> print(a)
a = 4   # a now refers to the new object created by the function

Passing mutable objects:
def increment(n):
    n.append([4])

>>> L = [1, 2, 3]
>>> increment(L)
>>> print(L)
L = [1, 2, 3, 4]   # L changed!

Variable 'L' refers to the list object containing references/pointers to the
three immutable integers.
Inside the function, the variable 'n' is first assigned to the same list 
object as 'L' but the 'append' method modifies the list in place and alters
the list object that both variables refer to.

Another example:
def assign_value(n, v):
    n = v

>>> L1 = [1, 2, 3]
>>> L2 = [4, 5, 6]
>>> assign_value(L1, L2)
>>> print(L1)
[1, 2, 3]

The variable 'L1' did not change because the function acted upon the 'n' and 'v'
variables that were assigned to the same list objects as 'L1' and 'L2' inside
the function. The 'n' variable was reassigned from the list object that 'L1'
refers to to the list object that 'v' and 'L2' refer to.

def my_func(x, y)
    return x+y

print(my_func(8, 9))

'x' and 'y' refer to the integer objects 8 and 9, respectively. The names 'x' 
and 'y' are local to the function and are discarded after the function returns.
The values that those variables/names refered to might continue to exist if
other variables refer to them.
"""


import sys
from collections import Mapping, Container


def deep_getsizeof(o, ids):
    """
    Finds the memory footprint of a Python object.

    This is a recursive function that drills down a Python object graph
    like a dictionary holding nested dictionaries with lists of lists
    and tuples and sets.
    
    The sys.getsizeof function does a shallow size of only. It counts each
    object inside a container as pointer only regardless of how big it
    really is.
    
    :param o: the object
    :param ids: empty set to add objects ids
    :return: memory footprint of object in bytes
    """
    
    d = deep_getsizeof
    if id(o) in ids:
        return 0
    
    r = sys.getsizeof(o)
    ids.add(id(o))
    
    if isinstance(o, str) or isinstance(0, bytes):
        return r
    
    elif isinstance(o, Mapping):
        return r + sum(d(k, ids) + d(v, ids) for k, v in o.items())
    
    elif isinstance(o, Container):
        return r + sum(d(x, ids) for x in o)
    
    return r


def convert_bytes(o, ids):
    """
    Calls 'deep_getsizeof' function to calculate memory footprint of object in 
    bytes and returns string with size in bytes, Kilobytes, Megabytes and 
    Gigabytes.
    
    :param o: the object
    :param ids: empty set to add objects ids
    :return: memory footprint of object in bytes, Kilobytes, Megabytes and 
    Gigabytes.
    """
    
    r = deep_getsizeof(o, ids)
    
    bytes_size = r
    kilobytes_size = r / 1000
    megabytes_size = kilobytes_size / 1000
    gigabytes_size = megabytes_size / 1000
    
    return_string = str(bytes_size) + ' bytes \n' + \
                    str(kilobytes_size) + ' Kilobytes \n' + \
                    str(megabytes_size) + ' Megabytes \n' + \
                    str(gigabytes_size) + ' Gigabytes'

    return return_string


########
# Test #
########

#obj_list = [2,2.4,'12',[],set(),dict()]
#
#for o in obj_list:
#    print(convert_bytes(o,set()))
#    print('\n')
#
#print(convert_bytes(obj_list,set()))
    
