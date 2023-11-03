import Link from "next/link";
import { ItemContainer, ItemLabel } from "../styles";

function DetailItem(props: {
  label: string;
  children?: string;
  link?: string;
  onClick?: () => void;
}) {
  return props.children ? (
    <ItemContainer>
      <ItemLabel>{props.label}</ItemLabel>
      {props.link ? (
        <Link href={props.link} passHref onClick={props.onClick}>
          {props.children}
        </Link>
      ) : (
        props.children
      )}
    </ItemContainer>
  ) : null;
}

export default DetailItem;
