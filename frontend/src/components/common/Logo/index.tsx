import Image from "next/image";
import { FC } from "react";
import logo from "src/common/images/explore-white.svg";

export const Logo: FC = () => {
  return (
    <Image
      alt="cellxgene logo"
      data-test-id="logo"
      height={24}
      src={logo}
      width={24}
    />
  );
};
