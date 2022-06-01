import Image from "next/image";
import { FC } from "react";
import logo from "src/common/images/explore-white.svg";

interface LogoProps {
  dataTestId: string;
}

export const Logo: FC<LogoProps> = (props) => {
  return (
    <Image
      alt="cellxgene logo"
      data-test-id={props.dataTestId}
      height={24}
      src={logo}
      width={24}
    />
  );
};
