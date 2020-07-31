import { SystemStyleObject } from "@styled-system/css";
import React, { FC } from "react";
import { Button, Text } from "theme-ui";

interface Props {
  label: string;
  disabled?: boolean;
  handleClick: () => void;
}

const StyledButton: FC<Props> = ({ label, disabled, handleClick }) => {
  return (
    <Button
      sx={
        {
          backgroundColor: disabled ? "muted" : "steelblue",
          borderRadius: "1px",
          color: disabled ? "gray" : "white",
          cursor: disabled ? "default" : "pointer",
          fontWeight: "100",
          height: "35px",
          margin: "0px 5px",
          minWidth: "100px",
          padding: "5px 20px",
          pointerEvents: disabled ? "none" : "default",
          textAlign: "center",
        } as SystemStyleObject
      }
      onClick={disabled ? undefined : handleClick}
    >
      <Text sx={{ fontSize: [0] }}>{label}</Text>
    </Button>
  );
};

export default StyledButton;
