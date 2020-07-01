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
          pointerEvents: disabled ? "none" : "default",
        } as SystemStyleObject
      }
      onClick={disabled ? undefined : handleClick}
    >
      <Text sx={{ fontSize: [0] }}>{label}</Text>
    </Button>
  );
};

export default StyledButton;
