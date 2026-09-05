#!/bin/bash
#
# give a hint about detrmining prefix length
#
echo "         1         2         3         4         5         6         7"
echo "1234567890123456789012345678901234567890123456789012345678901234567890"
sed '/^#/d;/^$/d;/^\$/d;' $1 | head -5
