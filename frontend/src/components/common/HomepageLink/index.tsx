import Link from "next/link";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Logo } from "../Logo";

export const HomepageLink: FC = () => {
  return (
    <Link href={ROUTES.HOMEPAGE} passHref>
      <a href="passHref">
        <Logo />
      </a>
    </Link>
  );
};
