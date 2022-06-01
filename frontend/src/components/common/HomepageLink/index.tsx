import Link from "next/link";
import { FC } from "react";
import { ROUTES } from "src/common/constants/routes";
import { Logo } from "../Logo";

interface HomepageLinkProps {
  dataTestId: string;
}

export const HomepageLink: FC<HomepageLinkProps> = (props) => {
  return (
    <Link href={ROUTES.HOMEPAGE} passHref>
      <a href="passHref">
        <Logo dataTestId={props.dataTestId} />
      </a>
    </Link>
  );
};
