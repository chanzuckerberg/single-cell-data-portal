import Image from "next/image";
import { FC } from "react";
import logo from "src/common/images/cellxgene-discover-logo.svg";

export const Logo: FC = () => {
  return (
    <Image
      alt="CELLxGENE logo"
      data-testid="logo"
      height={28}
      src={logo}
      width={28}
    />
  );
};
