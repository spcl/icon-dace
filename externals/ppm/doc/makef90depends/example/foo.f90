MODULE foo
  TYPE rhubarb
    INTEGER :: i
  END TYPE rhubarb
END MODULE foo

MODULE banana
  USE foo, ONLY: rhubarb
END MODULE banana
