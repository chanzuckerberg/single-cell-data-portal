import Link from "next/link";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Logo } from "../Logo";

interface Props {
  homeUrl?: string;
}

export const HomepageLink: FC<Props> = ({ homeUrl }) => {
  return (
    <Link href={homeUrl || ROUTES.HOMEPAGE} passHref>
      <a href="passHref">
        <Logo />
      </a>
    </Link>
  );
};
