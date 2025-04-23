import Head from "next/head";
import Image from "next/image";
import { ROUTES } from "src/common/constants/routes";
import rawCellxgeneLogo from "src/components/common/staticPages/cellxgene.png";
import {
  CommonStyle,
  Layout,
  TOSStyle,
} from "src/components/common/staticPages/style";

const ToS = (): JSX.Element => {
  return (
    <Layout>
      <CommonStyle>
        <TOSStyle>
          <Head>
            <title>Terms of Service - CZ CELLxGENE</title>
          </Head>
          <header>
            <Image
              data-testid="cellxgene-logo"
              src={rawCellxgeneLogo}
              alt="CELLxGENE logo"
              width="119"
              height="35"
            />
          </header>

          <br />

          <h1 id="tos">Terms of Use</h1>

          <p>Last Updated: April 23, 2025.</p>

          <p>
            <span className="caps">
              ANY DISPUTE BETWEEN YOU AND US IS SUBJECT TO BINDING ARBITRATION,
              AS SET FORTH IN <a href="#section6">SECTION 6</a>. PLEASE READ THE
              ARBITRATION PROVISION IN THIS AGREEMENT AS IT AFFECTS YOUR RIGHTS
              UNDER THIS CONTRACT.
            </span>
          </p>

          <br />

          <h2>Introduction</h2>

          <p>
            These Terms of Use (the “<strong>Terms</strong>”) apply to everyone
            who accesses and uses the version of CELLxGENE located at
            cellxgene.cziscience.com (the “<strong>Site</strong>”; those
            individuals who access and use it, “<strong>you</strong>”).
          </p>

          <p>
            The purpose of the Site is to provide a version of the interactive
            data explorer called CELLxGENE that enables fast visualizations of
            curated single-cell transcriptomics datasets (those datasets, the “
            <strong>Data</strong>”). By leveraging modern web development
            techniques to enable those visualizations, we hope to enable
            biologists and computational researchers to explore that Data.
          </p>

          <p>
            You agree that by accessing the Site, you are entering into a
            legally-binding contract with the Chan Zuckerberg Initiative
            Foundation, a 501(c)(3) nonprofit private foundation (
            <strong>“Provider,” “we,” “us”</strong>) and agree to be bound by
            these Terms as well as our{" "}
            <a href={ROUTES.PRIVACY}>Privacy Policy</a>. If you don’t agree to
            these Terms, or to our <a href={ROUTES.PRIVACY}>Privacy Policy</a>,
            don’t access the Site.
          </p>

          <p>
            Please read our full Terms and{" "}
            <a href={ROUTES.PRIVACY}>Privacy Policy</a> for complete details,
            but here is the key information you should know:
          </p>

          <ul>
            <li>
              This is a <strong>free service</strong> we provide in order to
              advance biomedical research. The Site should not be relied upon to
              provide medical advice, diagnosis, or treatment.
            </li>
            <li>
              The data made available in the service is{" "}
              <strong>publicly available</strong> and is{" "}
              <strong>not personally identifiable</strong>.
            </li>
            <li>
              We use the privacy-friendly Plausible service to collect{" "}
              <strong>basic analytics</strong> about Site traffic (e.g. the
              number of visitors and page views) so we know how it’s being used.
              This helps us improve the service as well as understand how it’s
              used.
            </li>
          </ul>

          <p className="caps">
            PLEASE READ THESE TERMS CAREFULLY AS USE OF THE SERVICE CONSTITUTES
            ACCEPTANCE OF THESE TERMS AND CONDITIONS.
          </p>

          <ol>
            <li>
              <h4>Data Submissions.</h4> If you would like to submit Data to the
              Site, please contact{" "}
              <a href="mailto:cellxgene@chanzuckerberg.com">
                cellxgene@chanzuckerberg.com
              </a>
              . Keep in mind that – if the Data you wish to submit is not
              already published in a public archive such as GEO, as an open link
              in a publication, or in a GitHub repository – you will be required
              to grant CZIF permission to use, display and create derivative
              works (e.g. visualizations) of the Data for purposes of offering
              the Site, and must therefore have the authority to give that
              permission. You can find the data submission policy and
              documentation for contributing data{" "}
              <a
                href="https://cellxgene.cziscience.com/docs/032__Contribute%20and%20Publish%20Data"
                target="_blank"
              >
                here
              </a>
              .
            </li>
            <li>
              <h4>Site Limitations.</h4> We may restrict or terminate your
              access to the Site at any time, without notice, and for any reason
              including for breach of these Terms. The Data on the Site has been
              compiled from a variety of sources, and is subject to change
              without notice. CZIF does not investigate, monitor or check Data
              for accuracy, appropriateness, completeness, or other reliability.
              As such, you agree that CZIF shall not be responsible for any Data
              or failure to include any Data or updates thereto. Your use of the
              Data is at your own risk.
            </li>
            <li>
              <h4>Our Rights in the Services.</h4> We reserve all rights, title,
              and interest in the Site. Subject to your compliance with these
              Terms, we grant you a limited right to access and use the Site.
              Using the Site does not give you any right, title, or interest in
              our Site, other than the right we explicitly grant you herein. The
              trademarks, service marks, graphics, and logos used for our
              Services, whether registered or unregistered, are owned by us or
              our licensors.
            </li>
            <li>
              <h4>Service Availability.</h4>We reserve the right to modify,
              update, pause, suspend, or discontinue the Site at any time. We
              strive to maintain continuous availability but make no guarantees
              regarding the availability of the Site.
            </li>
            <li>
              <h4>Feedback.</h4>We welcome feedback to improve the Site. By
              providing feedback, you assign to us all rights, title, and
              interest in the feedback, with no entitlement to compensation or
              rights to any resulting improvements. If such assignment is not
              permitted by law, you grant us a non-exclusive, irrevocable,
              transferable, sub-licensable, royalty-free, worldwide, perpetual
              license to use the feedback.
            </li>
            <li>
              <h4>Monitoring and Removal.</h4>We reserve the right to monitor
              and audit your use of the Site, including your submissions, to
              ensure compliance with these Terms. We may remove or restrict
              access to any of your submissions at our discretion if they
              violate these Terms or are otherwise objectionable.
            </li>
            <li>
              <h4>Third Party Materials</h4>We do not warrant or endorse and do
              not assume and will not have any liability or responsibility to
              you or any other person for any third-party materials we make
              available to you via the Site, or for any other materials,
              products, or services of third parties. Third-party materials or
              services are provided solely as a convenience to you.
            </li>
            <li>
              <h4>Copyright Complaints and DMCA.</h4>We respect the intellectual
              property rights of others. We require that content posted by you
              does not violate the intellectual property rights of third
              parties. Please see our{" "}
              <a href="/dmca/">
                Digital Millennium Copyright Act (DMCA) Policy
              </a>{" "}
              for more information. If you believe your intellectual property
              rights have been violated through the Services, information on how
              to contact us is available in our DMCA Policy.
            </li>
            <li>
              <h4>Disclaimer.</h4>{" "}
              <span className="caps">
                THE SITE AND THE DATA ARE PROVIDED “AS IS” AND “AS AVAILABLE”
                WITH ALL FAULTS, AND WE, OUR PARENTS, AFFILIATES, RELATED
                COMPANIES, OFFICERS, DIRECTORS, EMPLOYEES, AGENTS,
                REPRESENTATIVES, PARTNERS, AND LICENSORS (THE “PROTECTED
                PARTIES”) HEREBY DISCLAIM ALL WARRANTIES AND CONDITIONS, WHETHER
                EXPRESS OR IMPLIED, INCLUDING WARRANTIES OF MERCHANTABILITY,
                FITNESS FOR A PARTICULAR PURPOSE, OR NON-INFRINGEMENT. THE
                PROTECTED PARTIES MAKE NO WARRANTY OR REPRESENTATION AND
                DISCLAIM ALL RESPONSIBILITY AND LIABILITY FOR: (A) THE
                COMPLETENESS, ACCURACY, AVAILABILITY, TIMELINESS, SECURITY, OR
                RELIABILITY OF THE SITE; (B) ANY HARM TO YOUR COMPUTER SYSTEM,
                LOSS OF DATA, OR OTHER HARM THAT RESULTS FROM YOUR ACCESS TO OR
                USE OF THE SITE; (C) THE OPERATION OR COMPATIBILITY WITH ANY
                OTHER APPLICATION OR PARTICULAR SYSTEM OR DEVICE; (D) WHETHER
                THE SITE WILL MEET YOUR REQUIREMENTS OR BE AVAILABLE ON AN
                UNINTERRUPTED, SECURE, OR ERROR-FREE BASIS; AND (E) THE DELETION
                OF, OR FAILURE TO STORE OR TRANSMIT, YOUR SUBMISSIONS, AND OTHER
                COMMUNICATIONS MAINTAINED BY THE SITE. NO ADVICE OR INFORMATION,
                WHETHER ORAL OR WRITTEN, OBTAINED FROM THE PROTECTED PARTIES OR
                THROUGH THE SITE WILL CREATE ANY WARRANTY OR REPRESENTATION NOT
                EXPRESSLY MADE HEREIN.
              </span>
            </li>

            <li>
              <h4>Limitation of Liability.</h4>{" "}
              <span className="caps">
                TO THE MAXIMUM EXTENT PERMITTED BY APPLICABLE LAW, THE PROTECTED
                PARTIES WILL NOT BE LIABLE FOR ANY INDIRECT, INCIDENTAL,
                SPECIAL, PUNITIVE, OR CONSEQUENTIAL DAMAGES OF ANY KIND
                (INCLUDING LOST PROFITS, LOST DATA, BUSINESS INTERRUPTION, OR
                LOSS OF GOODWILL) IRRESPECTIVE OF WHETHER SUCH DAMAGES ARISE
                FROM CLAIMS BROUGHT IN CONTRACT, TORT, NEGLIGENCE, WARRANTY,
                STRICT LIABILITY, OR ANY OTHER THEORY AT LAW OR IN EQUITY, AND
                EVEN IF ANY PROTECTED PARTY HAS BEEN ADVISED OF THE POSSIBILITY
                OF SUCH DAMAGES. WITHOUT LIMITING THE FOREGOING, TO THE MAXIMUM
                EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT WILL THE
                PROTECTED PARTIES’ AGGREGATE LIABILITY ARISING OUT OF OR
                RELATING TO THESE TERMS OR THE SERVICE EXCEED USD $100. THE
                FOREGOING LIMITATIONS WILL APPLY EVEN IF THE ABOVE STATED REMEDY
                FAILS OF ITS ESSENTIAL PURPOSE.
              </span>
            </li>

            <li>
              <h4>Indemnification.</h4>
              <ol className="section5">
                <li>
                  You shall indemnify, defend and hold the Protected Parties
                  harmless from and against, and shall pay all damages, costs,
                  fees and expenses (including reasonable attorneys’ fees and
                  expenses) relating to, any third party (including government
                  entity) claim, action, suit or other proceeding (a “Claim”) to
                  the extent arising from: (1) your gross negligence, willful
                  misconduct or fraud; (2) your data submissions; (3) your
                  violation or breach of these Terms; (4) your violation of any
                  rights of any third party; and/or (5) any misrepresentation
                  you make regarding your permission to submit data to CELLxGENE
                  for public use.
                </li>
                <li>
                  Section 12.1 indemnification is conditioned upon, CZIF giving
                  you written notice of any such Claim, and giving you control
                  of the defense and settlement of any such Claim, and
                  cooperating with you in such defense. Notwithstanding anything
                  to the contrary, (1) CZIF may participate in defense of such
                  Claim with its own counsel at its own expense and (2) you may
                  not settle any Claim without CZIF’s prior written consent,
                  which will not be unreasonably withheld, unless it
                  unconditionally releases CZIF of all liability, obligation,
                  and fault.
                </li>
              </ol>
            </li>

            <li id="section6">
              <h4>Arbitration Agreement and Class Action Waiver.</h4>
              <ol className="section6">
                <li>
                  <h5>Applicability.</h5> In the unlikely event we end up in a
                  legal dispute, you and CZIF agree that all Disputes, including
                  Enforceability Disputes, will be resolved exclusively in
                  binding arbitration on an individual basis, except that you
                  and CZIF are not required to arbitrate IP Disputes.
                  Notwithstanding the foregoing, either you or CZIF may bring an
                  individual action in small claims court.
                  <ol>
                    <li>
                      (1) A <strong>“Dispute”</strong> means a dispute, claim or
                      controversy arising out of or relating to the Site or
                      these Terms; whether that dispute is (1) based on past,
                      present or future events; or (2) in tort, contract,
                      warranty, state, regulation, or other legal or equitable
                      basis.
                    </li>
                    <li>
                      (2) An
                      <strong>“Enforceability Dispute”</strong> means a Dispute
                      relating to the interpretation, applicability, or
                      enforceability of this Arbitration Agreement, including
                      the formation of the contract, the arbitrability of any
                      Dispute, and any claim that all or any part of this
                      Arbitration Agreement is void or voidable.
                    </li>
                    <li>
                      (3) An <strong>“IP Dispute”</strong> means a Dispute
                      relating to the ownership or enforcement of intellectual
                      property rights.
                    </li>
                  </ol>
                </li>
                <li>
                  <h5>Waivers.</h5>
                  <ol>
                    <li>
                      <strong>Waiver of Jury Right.</strong> YOU AND CZIF ARE
                      EXPRESSLY GIVING UP ALL RIGHTS TO A JURY TRIAL OR COURT
                      TRIAL BEFORE A JUDGE, EXCEPT AS EXPRESSLY PROVIDED IN THIS
                      ARBITRATION AGREEMENT. The arbitrator’s decision will be
                      final and binding on both you and us, subject to review
                      solely on the grounds set forth in the Federal Arbitration
                      Act (“FAA”).
                    </li>
                    <li>
                      <strong>Waiver of Class or Consolidated Actions.</strong>
                      YOU AND CZIF AGREE THAT ALL DISPUTES MUST BE ARBITRATED OR
                      LITIGATED ON AN INDIVIDUAL BASIS AND NOT ON A CLASS,
                      COLLECTIVE ACTION, OR REPRESENTATIVE BASIS. The validity
                      of this waiver – and whether an action may proceed as a
                      class, collective, or representative action – must be
                      decided by a court.
                    </li>
                  </ol>
                </li>
                <li>
                  <h5>Initiating a Dispute.</h5> To initiate a Dispute, a party
                  must send to the other party written notice of that Dispute
                  containing: (a) the name, address, and contact information of
                  the party giving notice; (b) the facts giving rise to the
                  Dispute; and (c) the relief requested. Notices sent to CZIF
                  must be sent by mail to the address provided in Section 18
                  below.{" "}
                  <p>
                    You and we agree that the parties shall (in good faith) meet
                    and attempt to resolve the Dispute within 30 days. If the
                    Dispute is not resolved during that time period, then you
                    and a representative of CZIF shall (in good faith) meet and
                    attempt to resolve the Dispute through non-binding mediation
                    with a mutually agreed-upon mediator within 30 additional
                    days. If you and we do not reach an agreement to resolve the
                    dispute within that 60-day period, you or we may commence an
                    arbitration proceeding or file a claim in small claims
                    court.
                  </p>
                </li>
                <li>
                  <h5>Arbitration Rules and Procedure.</h5>
                  <ol>
                    <li>
                      <strong>Rules.</strong> The Federal Arbitration Act
                      governs the interpretation and enforcement of this
                      Arbitration Agreement. Judicial Arbitration & Mediation
                      Services, Inc. (“JAMS”) will administer the arbitration
                      before a single arbitrator, and the arbitration will be
                      initiated and conducted according to the Streamlined
                      Arbitration Rules and Procedures (the “JAMS Rules”), to
                      the extent they are not inconsistent with the terms of
                      this Arbitration Agreement. The JAMS Rules and
                      instructions about how to initiate an arbitration are
                      available at
                      https://www.jamsadr.com/rules-streamlined-arbitration (as
                      of the date of this agreement) or 1-800-352-5267.
                    </li>
                    <li>
                      <strong>Fees.</strong> Pursuant to the JAMS Consumer
                      Arbitration Minimum Standards, CZIF will bear all costs of
                      the arbitration (including any JAMS Case Management Fee
                      and all professional fees for the arbitrator’s services),
                      except for the filing fee if you are the party initiating
                      the arbitration.
                    </li>
                    <li>
                      <strong>Manner and Location of Arbitration.</strong> You
                      may choose to have the arbitration conducted by telephone,
                      in writing, online, or in person. If in person, you may
                      choose to have the arbitration conducted (a) in San Mateo
                      County, California, (b) if you are not domiciled in the
                      United States, in your country of domicile, or (c) at
                      another location that you and we agree upon.
                    </li>
                  </ol>
                </li>
                <li>
                  <h5>Opt out.</h5> You may opt out of this Arbitration
                  Agreement by notifying us no later than 30 days after first
                  becoming subject to it. Your notice must include your name,
                  address, and a clear statement that you want to opt out of
                  this Arbitration Agreement. Notices sent to CZIF must be sent
                  by mail to the address provided in Section 18 of this
                  Agreement.
                </li>
                <li>
                  <h5>Severability.</h5> If any portion of this Arbitration
                  Agreement is found to be unlawful, void or for any reason
                  unenforceable, then that portion shall be severed and the
                  remainder of this Arbitration Agreement shall be given full
                  force and effect.
                </li>
              </ol>
            </li>

            <li id="section7">
              <h4>General Terms.</h4>
              <ol className="section7">
                <li>
                  These Terms constitute the entire agreement between the
                  parties with respect to the subject matter hereof and your use
                  of the Site, and supersedes all other agreements and
                  understandings, both written and oral, between the parties
                  with respect to the subject matter hereof.
                </li>
                <li>
                  If any provision in these Terms is held invalid or
                  unenforceable, then that provision shall be deemed severable
                  from these Terms and shall not affect the validity and
                  enforceability of any remaining provisions.
                </li>
                <li>
                  No waiver by either party of any breach or default hereunder
                  shall be deemed to be a waiver of any preceding or subsequent
                  breach or default.
                </li>
                <li>
                  You may not assign or transfer any rights or obligations under
                  these Terms without our consent. However, you agree that we
                  may assign these Terms in connection with a reorganization, or
                  to a successor or assign that agrees to assume our obligations
                  under these Terms (and{" "}
                  <a href={ROUTES.PRIVACY}>Privacy Policy</a>) without your
                  consent.
                </li>
              </ol>
            </li>
            <li id="section8">
              <h4>Governing Law.</h4> These Terms and any Dispute between you
              and CZIF will be governed by California law and/or applicable
              federal law (including the Federal Arbitration Act) without regard
              to its choice of law or conflicts of law principles.
            </li>
            <li id="section9">
              <h4>Jurisdiction and Venue.</h4> Subject to and without waiver of
              the arbitration provisions in Section 13, you agree that any
              judicial proceedings (other than small claims actions) will be
              brought in and you hereby consent to the exclusive jurisdiction
              and venue in the state courts in the city and county of San Mateo,
              California, or federal court for the Northern District of
              California. For countries where this is not permissible, this
              won’t deprive you of any protection you have under the law of the
              country where you live, or access to the courts in that country.
            </li>
            <li id="section10">
              <h4>Update to Terms.</h4> We may modify these Terms from time to
              time in which case we will update the “Last Updated” date at the
              top of these Terms. If we make changes that are material, we will
              use reasonable efforts to attempt to notify you, such as by e-mail
              and/or by placing a notice on the Site. However, it is your sole
              responsibility to review these Terms from time to time to view any
              such changes. The updated Terms will be effective as of the time
              of posting, or such later date as may be specified in the updated
              Terms. Your continued access or use of the Services after the
              modifications have become effective will be deemed your acceptance
              of the modified Terms. No amendment shall apply to a Dispute for
              which an arbitration has been initiated prior to the change in
              Terms.
            </li>
            <li>
              <h4>How To Contact Us.</h4> Notice under these Terms must be in
              writing and deemed to have been given on the date delivered by a
              nationally recognized express mail service, such as Federal
              Express, or by certified and registered mail (signature for
              receipt required) to CZIF as follows:
              <br />
              <br />
              <address>
                Chan Zuckerberg Initiative <br />
                2682 Middlefield Road, Suite i <br />
                Redwood City, CA 94063
                <br />
                Attn: General Counsel <br />
                Email: courtesy copy:{" "}
                <a href="mailto:legalczi1@chanzuckerberg.com">
                  legalczi1@chanzuckerberg.com
                </a>{" "}
                (email does not constitute notice)
              </address>
            </li>
          </ol>
        </TOSStyle>
      </CommonStyle>
    </Layout>
  );
};

export default ToS;
