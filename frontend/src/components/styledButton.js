import React from "react"
import { Button, Text } from "theme-ui"

const StyledButton = ({ label, disabled, onclick }) => {
  return (
    <Button sx={{
                  height: "35px",
                  minWidth: "100px",
                  padding: "5px 20px",
                  textAlign: "center",
                  borderRadius: "1px",
                  fontWeight: "100",
                  margin: "0px 5px",
                  color: disabled ? "gray" : "white",
                  backgroundColor: disabled ? "muted" : "steelblue",
                  cursor: disabled ? "default" : "pointer",
                  pointerEvents: disabled ? "none" : "default"
                }}
             onClick={disabled ? undefined : onclick}>
      <Text sx={{ fontSize: [0] }}>{label}</Text>
    </Button>
  )
}

export default StyledButton
