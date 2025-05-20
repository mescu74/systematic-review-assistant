# Process lines that consist only of a plain code fence (```)
# possibly surrounded by whitespace.
/^[[:space:]]*\`\`\`[[:space:]]*$/{
  # Exchange pattern space (current line) with hold space (state flag).
  # PS now holds the flag, HS holds the current line (the fence).
  x

  # Check if the flag (now in PS) is empty.
  # An empty flag means we were OUTSIDE a code block.
  # So, the current fence line in HS is an OPENING fence.
  /^$/{
    # Set the flag to "M" (marker for "inside a block").
    s/.*/M/
    # Swap back: PS gets the fence line, HS gets the new flag "M".
    x
    # Modify the fence line in PS: add 'python'.
    # This regex captures leading whitespace (\1), the backticks (\2),
    # and trailing whitespace (\3), then reassembles them with 'python'.
    s~^\([[:space:]]*\)\(\`\`\`\)\([[:space:]]*\)$~\1\2python\3~
    # Print the modified opening fence.
    p
    # Branch to the end of the script for this line (skip default action).
    b
  }

  # If the flag (in PS) was not empty, it means it was "M" (or our marker).
  # This indicates we were INSIDE a code block.
  # So, the current fence line in HS is a CLOSING fence.
  # PS currently holds "M".
  # Clear the flag (set PS to empty, marking "outside a block").
  s/M//
  # Swap back: PS gets the fence line, HS gets the new empty flag.
  x
  # Print the original closing fence (unmodified).
  p
  # Branch to the end of the script for this line.
  b
}

# Default action for any line that was NOT a plain ` ``` ` fence
# (e.g., text, content inside a code block, or ```javascript).
# These lines are printed as is. The hold space (state flag) remains unchanged.
p
