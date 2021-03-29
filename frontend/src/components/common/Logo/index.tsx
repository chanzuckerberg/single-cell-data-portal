import Image from "next/image";
import React, { FC } from "react";
import logo from "src/common/images/logo.svg";
import { LogoWrapper } from "./style";

export const Logo: FC = () => {
  return (
    <LogoWrapper data-test-id="logo">
      <Image width="147px" height="35px" src={String(logo)} />
    </LogoWrapper>
  );
};
