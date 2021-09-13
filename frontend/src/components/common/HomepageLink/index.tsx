import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Logo } from "../Logo";
import { Link } from "./style";

export const HomepageLink: FC = () => {
  return (
    <Link href={ROUTES.HOMEPAGE} passHref>
      <a href="passHref">
        <Logo />
      </a>
    </Link>
  );
};
